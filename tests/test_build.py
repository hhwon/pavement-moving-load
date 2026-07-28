"""
Offline checks for everything up to "the run directory is ready".

No Abaqus, no cluster, no network — deck generation is pure text, so it can
be verified anywhere:

    python tests/test_build.py
    # or:  python -m pytest tests/ -q
"""
from __future__ import annotations

import json
import re
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from movingload.config import (MOTIONS, STRUCTURES, CaseConfig, Mesh, Motion,
                               Step, Structure, Tyre, make_case)
from movingload.dload import source as dload_source
from movingload.inp_builder import build, deck_lines, metadata
from movingload.mesh import build_grid
from movingload.slurm import Resources, local_script, slurm_script
from movingload import submit as submit_mod


def _tmp() -> Path:
    return Path(tempfile.mkdtemp())


# --------------------------------------------------------------- config

def test_presets_are_valid():
    """Every shipped preset must survive its own validation."""
    for name, s in STRUCTURES.items():
        assert s.name == name
        thick = s.thicknesses(Mesh().LY)
        assert len(thick) == len(s.E) == len(s.nu) == len(s.names) == len(s.n_elems)
        assert thick[-1] > 0, f"{name}: subgrade has no room left"
    for name, m in MOTIONS.items():
        assert m.t_final > 0
        assert (name == "vertical") == (not m.horizontal)


def test_structure_rejects_ragged_lists():
    try:
        Structure(name="bad", description="", h=[0.1, 0.2],
                  E=[1e9, 1e9], nu=[0.3, 0.3],       # one short: no subgrade
                  names=["A", "B"], n_elems=[2, 2])
    except ValueError as exc:
        assert "expected 3" in str(exc)
        return
    raise AssertionError("ragged layer lists were accepted")


def test_structure_rejects_overfull_domain():
    s = Structure(name="deep", description="", h=[3.0, 3.0],
                  E=[1e9] * 3, nu=[0.3] * 3, names=list("ABC"), n_elems=[2, 2, 2])
    try:
        s.thicknesses(5.0)          # 6 m of finite layers in a 5 m domain
    except ValueError as exc:
        assert "do not fit" in str(exc)
        return
    raise AssertionError("layers deeper than the domain were accepted")


def test_tyre_must_clear_the_symmetry_plane():
    try:
        Tyre(x_centre=0.05, radius=0.1065)
    except ValueError as exc:
        assert "symmetry plane" in str(exc)
        return
    raise AssertionError("a tyre print crossing x=0 was accepted")


def test_tyre_force_and_spacing():
    t = Tyre()
    assert abs(t.force - 3.14159265358979 * t.radius ** 2 * t.pressure) < 1e-6
    assert abs(t.force / 1e3 - 24.94) < 0.01, "BZZ-100 quarter axle ~25 kN"
    assert abs(t.dual_spacing - 0.32) < 1e-12


def test_motion_timing_is_self_consistent():
    m = MOTIONS["vertical"]
    z_mid = 3.5
    t_star = m.time_over(z_mid)
    assert abs(m.z_at(t_star) - z_mid) < 1e-12
    assert abs(t_star - 0.5) < 1e-12, "mid-span crossing should be t = 0.5"


def test_motion_sign_validation():
    for bad in (dict(mu=-0.1), dict(mu_sign=0.0), dict(t_final=0.0)):
        try:
            Motion(**bad)
        except ValueError:
            continue
        raise AssertionError(f"Motion accepted {bad}")


def test_step_must_divide_the_motion():
    try:
        CaseConfig(case_id="x", structure=STRUCTURES["semi_rigid"],
                   motion=Motion(t_final=1.0), step=Step(dt=0.03))
    except ValueError as exc:
        assert "whole number" in str(exc)
        return
    raise AssertionError("a dt that does not divide t_final was accepted")


def test_mesh_must_resolve_the_tyre():
    try:
        CaseConfig(case_id="x", structure=STRUCTURES["semi_rigid"],
                   tyre=Tyre(x_centre=0.40, radius=0.1065))   # print to 0.507 m
    except ValueError as exc:
        assert "fine x band" in str(exc)
        return
    raise AssertionError("a tyre print outside the fine mesh band was accepted")


# ----------------------------------------------------------------- mesh

def test_grid_size_and_ids():
    g = build_grid(make_case("semi_rigid", "vertical"))
    assert (g.NX, g.NY, g.NZ) == (24, 31, 70)
    assert g.n_elements == 24 * 31 * 70 == 52_080
    assert g.node_id(0, 0, 0) == 1
    assert g.elem_id(0, 0, 0) == 1
    # i is fastest, then j, then k.
    assert g.elem_id(1, 0, 0) == 2
    assert g.elem_id(0, 1, 0) == 1 + g.NX
    assert g.elem_id(0, 0, 1) == 1 + g.NX * g.NY


def test_all_structures_share_the_grid_size():
    """Layer subdivisions were chosen so the three structures stay comparable."""
    sizes = {name: build_grid(make_case(name)).n_elements for name in STRUCTURES}
    assert len(set(sizes.values())) == 1, f"grids diverged: {sizes}"


def test_observation_column_tracks_the_tyre():
    cfg = make_case("semi_rigid", "vertical")
    g = build_grid(cfg)
    # The observation column must be the one containing the tyre centre.
    assert g.x[g.i_obs] <= cfg.tyre.x_centre <= g.x[g.i_obs + 1]
    # And it must sit at the bottom of the surface layer.
    assert g.j_obs == cfg.structure.n_elems[0] - 1
    assert abs(g.y[g.j_obs + 1] - cfg.structure.h[0]) < 1e-12
    assert len(g.obs_elem_ids) == g.NZ


def test_layer_rows_sum_to_ny():
    for name, s in STRUCTURES.items():
        g = build_grid(make_case(name))
        assert len(g.layer_of_row) == g.NY, name
        for ilayer, n in enumerate(s.n_elems):
            assert g.layer_of_row.count(ilayer) == n, f"{name} layer {ilayer}"


# ------------------------------------------------------------------ inp

def test_deck_structure():
    cfg = make_case("semi_rigid", "vertical")
    g = build_grid(cfg)
    text = "\n".join(deck_lines(cfg, g))

    assert "*ELEMENT, TYPE=C3D8R" in text
    assert text.count("*MATERIAL, NAME=") == len(cfg.structure.names)
    assert text.count("*ELASTIC") == len(cfg.structure.names)
    for s in ("SYM_X", "FAR_X", "FAR_Y", "FAR_Z0", "FAR_ZL",
              "TOP_ELEMS", "OBS_ELEMS"):
        assert s in text, f"missing set {s}"
    assert "*STATIC, DIRECT" in text
    assert "TOP_ELEMS, P3NU, 1.0" in text          # DLOAD-driven pressure
    assert "E33," in text                          # the requested history output
    assert text.rstrip().endswith("*END STEP")


def test_vertical_deck_has_no_dsload():
    cfg = make_case("semi_rigid", "vertical")
    text = "\n".join(deck_lines(cfg, build_grid(cfg)))
    assert "TRVEC3NU" not in text and "*DSLOAD" not in text


def test_traction_deck_has_dsload():
    for motion in ("accel", "brake"):
        cfg = make_case("semi_rigid", motion)
        text = "\n".join(deck_lines(cfg, build_grid(cfg)))
        assert "*DSLOAD, OP=NEW" in text, motion
        assert "TOP_ELEMS, TRVEC3NU, 1.0" in text, motion


def test_every_element_belongs_to_exactly_one_layer():
    cfg = make_case("rigid", "vertical")
    g = build_grid(cfg)
    text = "\n".join(deck_lines(cfg, g))
    seen: set[int] = set()
    for name in cfg.structure.names:
        block = re.search(rf"\*ELSET, ELSET={name}\n((?:[\d, ]+\n)+)", text)
        assert block, f"no elset for {name}"
        ids = {int(v) for v in block.group(1).replace("\n", ",").split(",") if v.strip()}
        assert not (ids & seen), f"{name} overlaps a previous layer"
        seen |= ids
    assert len(seen) == g.n_elements, "layer elsets do not cover every element"


def test_build_writes_both_files():
    d = _tmp()
    info = build(make_case("flexible", "accel"), d)
    assert Path(info["inp"]).is_file() and Path(info["mesh_json"]).is_file()
    meta = json.loads(Path(info["mesh_json"]).read_text())
    assert meta["case_id"] == "flexible_accel"
    assert meta["motion"]["mu"] == 0.30 and meta["motion"]["mu_sign"] == 1.0
    assert len(meta["observation"]["OBS_ELEM_IDS"]) == meta["NZ"]
    assert len(meta["X_NODES"]) == meta["NX"] + 1


def test_metadata_matches_the_deck():
    """The sidecar the extractor reads must describe the deck that was written."""
    cfg = make_case("rigid", "brake")
    g = build_grid(cfg)
    meta = metadata(cfg, g)
    text = "\n".join(deck_lines(cfg, g))
    block = re.search(r"\*ELSET, ELSET=OBS_ELEMS\n((?:[\d, ]+\n)+)", text)
    in_deck = [int(v) for v in block.group(1).replace("\n", ",").split(",") if v.strip()]
    assert in_deck == meta["observation"]["OBS_ELEM_IDS"]


# ---------------------------------------------------------------- dload

_NUM = re.compile(r"(?<![A-Za-z0-9_])([+-]?\d*\.?\d+)[dD]([+-]?\d+)")


def _statements(text: str) -> list[str]:
    """Executable Fortran only: comments dropped, literals and spacing normalised."""
    out = []
    for line in text.splitlines():
        s = line.strip()
        if not s or (s.startswith("C") and (len(s) == 1 or not s[1].isalnum())):
            continue
        s = re.sub(r"\s+", " ", s)
        s = _NUM.sub(lambda m: repr(float(m.group(1)) * 10 ** int(m.group(2))), s)
        s = re.sub(r"\s*([(),=])\s*", r"\1", s)
        out.append(s)
    return out


def test_vertical_dload_has_no_utracload():
    cfg = make_case("semi_rigid", "vertical")
    src = dload_source(cfg.tyre, cfg.motion, cfg.case_id)
    assert "SUBROUTINE DLOAD" in src
    assert "UTRACLOAD" not in src, "a vertical case must not emit a traction routine"


def test_traction_dload_emits_both_routines():
    for motion, expect_sign in (("accel", 1.0), ("brake", -1.0)):
        cfg = make_case("semi_rigid", motion)
        src = dload_source(cfg.tyre, cfg.motion, cfg.case_id)
        assert "SUBROUTINE DLOAD" in src and "SUBROUTINE UTRACLOAD" in src
        stmts = _statements(src)
        assert f"PARAMETER(MU_SIGN={expect_sign})" in stmts, motion
        assert "PARAMETER(MU=0.3)" in stmts, motion


def test_dload_parameters_follow_the_config():
    cfg = make_case("semi_rigid", "vertical")
    stmts = _statements(dload_source(cfg.tyre, cfg.motion))
    for name, value in (("P0", cfg.tyre.pressure), ("R", cfg.tyre.radius),
                        ("X_TYRE", cfg.tyre.x_centre),
                        ("ZSTART", cfg.motion.z_start),
                        ("ZEND", cfg.motion.z_end),
                        ("TFINAL", cfg.motion.t_final)):
        assert f"PARAMETER({name}={value!r})" in stmts, f"{name} not carried through"


def test_dload_is_fixed_form_fortran():
    """Columns 1-6 are reserved; continuation markers go in column 6."""
    src = dload_source(*[make_case("semi_rigid", "brake").tyre,
                         make_case("semi_rigid", "brake").motion])
    for line in src.splitlines():
        if not line or line[0] in "C*!":
            continue
        assert len(line) <= 72, f"line exceeds column 72: {line!r}"
        if line[:5].strip():
            raise AssertionError(f"statement text in columns 1-5: {line!r}")


# ---------------------------------------------------------------- slurm

def test_slurm_script_essentials():
    txt = slurm_script("job1", "job1.inp", resources=Resources(cpus=8, mem="16g"))
    assert txt.startswith("#!/bin/bash")
    # Without `interactive` the batch script exits and the job is torn down.
    assert "interactive" in txt
    assert "user=dload.f" in txt
    assert "--cpus-per-task=8" in txt and "--mem=16g" in txt
    assert "abaqus python extract_strain.py job1.odb" in txt


def test_partition_directive_is_optional():
    assert "--partition" not in slurm_script("j", "j.inp",
                                             resources=Resources(partition=""))
    assert "--partition=compute" in slurm_script(
        "j", "j.inp", resources=Resources(partition="compute"))


def test_module_load_is_optional():
    txt = slurm_script("j", "j.inp", module_load="")
    assert "abaqus assumed on PATH" in txt
    txt = slurm_script("j", "j.inp", module_load="module load abaqus/2023")
    assert "module load abaqus/2023" in txt


def test_local_script_runs_the_same_body():
    txt = local_script("j", "j.inp", cpus=2)
    assert "#SBATCH" not in txt
    assert "interactive" in txt and "extract_strain.py" in txt


# ------------------------------------------------------- run directory

def test_prepare_makes_a_self_contained_run_dir():
    d = _tmp()
    out = submit_mod.prepare("semi_rigid", "brake", out_root=d)
    run = out["run_dir"]
    for name in ("semi_rigid_brake.inp", "semi_rigid_brake_mesh.json",
                 "dload.f", "extract_strain.py", "submit.sh"):
        assert (run / name).is_file(), f"missing {name}"
    assert (run / "submit.sh").stat().st_mode & 0o111, "submit.sh not executable"
    # The extractor must find its sidecar by the ODB name alone.
    meta_name = f"{out['config'].case_id}_mesh.json"
    assert (run / meta_name).is_file()


def test_prepare_local_writes_run_sh_instead():
    d = _tmp()
    out = submit_mod.prepare("flexible", "vertical", out_root=d, local=True)
    assert (out["run_dir"] / "run.sh").is_file()
    assert not (out["run_dir"] / "submit.sh").exists()


def test_extractor_is_pure_ascii():
    """Abaqus' Python 2.7 rejects a source file containing non-ASCII bytes
    unless it declares an encoding (PEP 263).  Typographic dashes and quotes
    creep into docstrings very easily, and the failure only shows up on the
    solver host, after the solve has already been paid for.  Staying ASCII
    sidesteps the whole question.
    """
    path = REPO / "movingload" / "extract_strain.py"
    raw = path.read_bytes()
    if raw.startswith(b"#!") or b"coding" in raw[:100]:
        return          # an explicit encoding declaration also settles it
    for lineno, line in enumerate(raw.split(b"\n"), 1):
        try:
            line.decode("ascii")
        except UnicodeDecodeError:
            bad = next(b for b in line if b > 127)
            raise AssertionError(
                f"extract_strain.py:{lineno} contains non-ASCII byte "
                f"0x{bad:02x} and declares no encoding — Abaqus Python 2.7 "
                f"will refuse to parse it: {line[:70]!r}")


def test_extractor_stays_within_abaqus_python():
    """extract_strain.py runs under Abaqus' frozen interpreter (2.7 in most
    releases), so it must avoid f-strings, pathlib and annotations.

    Checks the parsed code, not the text: the module docstring legitimately
    mentions the very constructs it avoids.
    """
    import ast
    path = REPO / "movingload" / "extract_strain.py"
    tree = ast.parse(path.read_text(encoding="utf-8"))

    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported |= {a.name.split(".")[0] for a in node.names}
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.split(".")[0])
    assert "pathlib" not in imported, "pathlib is not in Abaqus Python 2.7"

    for node in ast.walk(tree):
        assert not isinstance(node, ast.JoinedStr), "f-string in extract_strain.py"
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            assert node.returns is None, f"return annotation on {node.name}()"
            args = node.args
            for a in (list(args.args) + list(args.kwonlyargs)
                      + list(getattr(args, "posonlyargs", []))):
                assert a.annotation is None, f"annotation on {node.name}({a.arg})"
        assert not isinstance(node, ast.AnnAssign), "annotated assignment found"


# ------------------------------------------------------------ sanitised

def test_no_private_paths_or_hosts_anywhere():
    """This repo is public: no personal paths, hosts or credentials.

    The patterns are assembled from fragments so that this file does not
    itself contain the literals it forbids — the scan covers the test file
    too.
    """
    parts = [
        "/" + "home" + r"/[a-z][a-z0-9_-]*/",   # personal home directories
        "/" + "mnt" + "/" + "public",           # site-specific mount point
        "@" + "local" + "host",                 # forwarded SSH gateway
        r"\b\d{1,3}(?:\.\d{1,3}){3}\b",         # bare IPv4 address
        r"[\w.+-]+@[\w-]+\.[a-z]{2,}",          # e-mail address
    ]
    forbidden = re.compile("|".join(parts), re.IGNORECASE)
    checked = 0
    for path in list(REPO.rglob("*.py")) + list(REPO.rglob("*.md")) \
            + list(REPO.rglob("*.toml")) + list(REPO.rglob("*.f")):
        if ".git" in path.parts or "runs" in path.parts:
            continue
        checked += 1
        for n, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            # Documented placeholders are fine.
            if "example.org" in line or "example.com" in line:
                continue
            hit = forbidden.search(line)
            assert not hit, (f"{path.relative_to(REPO)}:{n} leaks "
                             f"{hit.group(0)!r}: {line.strip()[:80]}")
    assert checked > 5, "sanitisation scan found almost nothing to check"


def test_clear_results_targets_only_output():
    """A re-used run directory must not hand back the previous run's CSVs.

    rsync runs without --delete, so `push` clears the old output first.  This
    issues `rm -f` on the cluster, so check exactly what it will and will not
    touch — and that unsafe paths are refused rather than expanded.
    """
    from movingload.cluster import Cluster

    cl = Cluster(host="user@example.org", remote_root="/scratch/x")
    sent = []
    cl.run = lambda cmd, check=True: sent.append(cmd)

    for path in ("/scratch/me/movingload/case", "~/runs/case"):
        sent.clear()
        assert cl.clear_results(path) is True, path
        assert len(sent) == 1
        cmd = sent[0]
        assert cmd.startswith("rm -f ")
        globs = set()
        for target in cmd[len("rm -f "):].split():
            base = target.rsplit("/", 1)[1]
            assert "*" in base and base != "*", target
            globs.add(base)
        # Inputs must survive.  Compare whole globs: '*.f' is a *substring* of
        # '*.fil' but not the same pattern, and only '*.fil' is deleted.
        for keep in ("*.inp", "*.f", "*.sh", "*.py", "*.json"):
            assert keep not in globs, f"cleanup would delete inputs: {keep}"
        # Concretely: the deck, the subroutine and the sidecar are untouched.
        import fnmatch
        for name in ("case.inp", "dload.f", "submit.sh", "run.sh",
                     "extract_strain.py", "case_mesh.json"):
            hit = [g for g in globs if fnmatch.fnmatch(name, g)]
            assert not hit, f"{name} would be deleted by {hit}"

    for unsafe in ("/", "~", "/tmp", "~/a; rm -rf /", "~/a/*/c", "~/`whoami`"):
        sent.clear()
        assert cl.clear_results(unsafe) is False, f"ran on {unsafe!r}"
        assert not sent


def test_cluster_refuses_without_configuration():
    """No baked-in SSH target: an unconfigured Cluster must fail loudly."""
    import os
    from movingload.cluster import Cluster, ClusterNotConfigured
    saved = {k: os.environ.pop(k, None)
             for k in ("MOVINGLOAD_SSH_HOST", "MOVINGLOAD_REMOTE_ROOT")}
    try:
        cl = Cluster()
        assert cl.host is None
        for call in (lambda: cl.check(), lambda: cl.remote_dir_for("x")):
            try:
                call()
            except ClusterNotConfigured:
                continue
            raise AssertionError("unconfigured Cluster did not refuse")
    finally:
        for k, v in saved.items():
            if v is not None:
                os.environ[k] = v


if __name__ == "__main__":
    passed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
            passed += 1
    print(f"\nAll {passed} checks passed.")

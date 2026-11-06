"""
Reduce an axisymmetric step-response ODB to the response table R(r, s).
Runs on the solver host, under Abaqus' own interpreter:

    abaqus python extract_axi.py <job>.odb

It reads `<job>_mesh.json` (written by `movingload.axi_inp`) for the
observation row and the step timing, and writes one CSV beside the ODB:

    <job>_axi_response.csv    s, element_id, r, E11, E22, E33, E12

s is the time since the load was switched on: 0 for the end of the INSTANT
(*STATIC) step, then the creep time through the CREEP_nn (*VISCO) steps.
E11 is radial, E22 axial (vertical), E33 hoop, E12 the engineering r-y
shear -- Abaqus' axisymmetric components, tension positive.  Turning this
into a moving-load history is done afterwards, in ordinary Python 3, by
`python -m movingload.superpose`.

Same constraints as extract_strain.py: Abaqus 2020 and earlier ship Python
2.7, so no f-strings, no pathlib, no annotations, ASCII only.
"""
from __future__ import print_function

import csv
import json
import os
import sys

try:
    from odbAccess import openOdb
except ImportError:
    sys.stderr.write(
        "ERROR: run this through Abaqus, e.g. `abaqus python extract_axi.py "
        "job.odb`, so that odbAccess is importable.\n")
    raise

COMPONENTS = ('E11', 'E22', 'E33', 'E12')


def find_meta(odb_path):
    base = os.path.splitext(os.path.basename(odb_path))[0]
    out_dir = os.path.dirname(os.path.abspath(odb_path)) or '.'
    path = os.path.join(out_dir, base + '_mesh.json')
    if not os.path.isfile(path):
        raise SystemExit('ERROR: no ' + path + ' -- it is written next to the '
                         '.inp by movingload.axi_inp; copy it beside the ODB.')
    return base, out_dir, path


def step_starts(meta):
    """Creep time at the start of each CREEP_nn step, in double precision.
    (step.totalTime in the ODB is single precision: summing it with the
    frame time misplaces s by ~1e-7 s, enough to break the ordering.)"""
    starts = {}
    t = 0.0
    for name, dt, period in meta['time']['creep_steps']:
        starts[str(name)] = t
        t += period
    return starts


def snap(s, known):
    """Replace s by the nearest expected output time; None if none is close."""
    best = min(known, key=lambda x: abs(x - s))
    if abs(best - s) <= 1e-4 * max(abs(best), 1e-6):
        return best
    return None


def collect(odb_path, meta):
    """Return rows (s, element_id, r, E11, E22, E33, E12), the number of
    frames that carried no strain for some observation element, and the
    number of frames whose time matched no expected output time."""
    obs = meta['observation']
    eids = obs['OBS_ELEM_IDS']
    r_of = {}
    for idx, eid in enumerate(eids):
        r_of[eid] = obs['r_centres'][idx]
    static_name = meta['time']['static_step']
    static_period = meta['time']['static_period']
    starts = step_starts(meta)
    known = meta['time']['output_times']

    odb = openOdb(odb_path, readOnly=True)
    rows = []
    incomplete = 0
    unmatched = 0
    try:
        for step_name in odb.steps.keys():
            step = odb.steps[step_name]
            for index, frame in enumerate(step.frames):
                if step_name == static_name:
                    # Only the loaded state; frame 0 is the unloaded start.
                    if abs(frame.frameValue - static_period) > 1e-6 * static_period:
                        continue
                    s = 0.0
                else:
                    if index == 0:
                        continue    # repeats the previous step's last frame
                    if step_name not in starts:
                        continue
                    s = snap(starts[step_name] + frame.frameValue, known)
                    if s is None:
                        unmatched += 1
                        continue
                if 'E' not in frame.fieldOutputs:
                    continue
                field = frame.fieldOutputs['E']
                labels = list(field.componentLabels)
                try:
                    cols = [labels.index(c) for c in COMPONENTS]
                except ValueError:
                    raise SystemExit('ERROR: strain components %s, expected %s'
                                     % (labels, COMPONENTS))
                got = {}
                for v in field.values:
                    if v.elementLabel in r_of and v.elementLabel not in got:
                        got[v.elementLabel] = [v.data[c] for c in cols]
                if len(got) < len(eids):
                    incomplete += 1
                    if not got:
                        continue
                for eid in eids:
                    if eid in got:
                        rows.append([s, eid, r_of[eid]] + got[eid])
    finally:
        odb.close()
    rows.sort(key=lambda row: (row[0], row[2]))
    return rows, incomplete, unmatched


def write_csv(path, header, rows):
    if sys.version_info[0] >= 3:
        handle = open(path, 'w', newline='')
    else:
        handle = open(path, 'wb')
    try:
        writer = csv.writer(handle)
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)
    finally:
        handle.close()
    print('wrote ' + path)


def main():
    if len(sys.argv) != 2:
        print('Usage: abaqus python extract_axi.py <job>.odb')
        return 1
    odb_path = sys.argv[1]
    if not os.path.isfile(odb_path):
        print('ERROR: no such ODB: ' + odb_path)
        return 1
    base, out_dir, meta_path = find_meta(odb_path)
    meta = json.load(open(meta_path))
    obs = meta['observation']
    print('case %s (%s): observation row at depth %.4f m, %d radii'
          % (meta.get('case_id', base), meta.get('material', '?'),
             obs['depth'], len(obs['OBS_ELEM_IDS'])))

    rows, incomplete, unmatched = collect(odb_path, meta)
    if unmatched:
        print('WARNING: %d frames fell on no expected output time and were '
              'dropped. Was the deck edited after the metadata was written?'
              % unmatched)
    if incomplete:
        print('WARNING: %d frames lacked strain on part of OBS_ROW. Check that '
              'the deck requests E on ELSET=OBS_ROW every increment.' % incomplete)
    if not rows:
        print('ERROR: no strain collected -- nothing to write.')
        return 2
    n_s = len(set(row[0] for row in rows))
    print('collected %d samples at %d of %d expected times'
          % (len(rows), n_s, len(meta['time']['output_times'])))
    write_csv(os.path.join(out_dir, base + '_axi_response.csv'),
              ['s', 'element_id', 'r'] + list(COMPONENTS), rows)
    return 0


if __name__ == '__main__':
    sys.exit(main())

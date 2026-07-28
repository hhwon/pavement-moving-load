"""
Reduce an Abaqus ODB to the longitudinal-strain history.  Runs on whatever
machine solved the case, under Abaqus' own interpreter:

    abaqus python extract_strain.py <job>.odb

It reads `<job>_mesh.json` (written next to the deck by
`movingload.inp_builder`) to learn which elements form the observation
column, then writes two CSVs beside the ODB:

    <job>_e33_history.csv    time, element_id, z_centre, E33   (whole column)
    <job>_e33_centre.csv     time, E33                         (mid-span only)

The first is the full (t, z) field on the observation column — enough to
draw the strain wave travelling with the tyre.  The second is the canonical
"approach – peak – leave" pulse at mid-span, which is what most analyses
actually compare.

Abaqus ships a frozen Python (2.7 in most releases, 3.x in recent ones), so
this file stays inside their common subset: no f-strings, no pathlib, no
type annotations, no dict comprehensions with complex targets.  Everything
else in this package targets Python 3.10+.
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
        "ERROR: run this through Abaqus, e.g. `abaqus python extract_strain.py "
        "job.odb`, so that odbAccess is importable.\n")
    raise


def find_meta(odb_path):
    """Locate the sidecar written next to the deck."""
    base = os.path.splitext(os.path.basename(odb_path))[0]
    out_dir = os.path.dirname(os.path.abspath(odb_path)) or '.'
    candidates = [
        os.path.join(out_dir, base + '_mesh.json'),
        os.path.join(out_dir, 'mesh.json'),
    ]
    for path in candidates:
        if os.path.isfile(path):
            return base, out_dir, path
    raise SystemExit(
        'ERROR: no mesh metadata found. Looked for:\n  ' +
        '\n  '.join(candidates) +
        '\nIt is written next to the .inp by movingload.inp_builder; copy it '
        'into the run directory alongside the ODB.')


def collect(odb_path, obs_eids, z_centres):
    """Pull E33 history for the observation elements.

    Returns a list of (time, element_id, z_centre, E33), and the number of
    observation elements that carried no E33 output (0 is the healthy case).
    """
    z_of_eid = {}
    for idx, eid in enumerate(obs_eids):
        z_of_eid[eid] = z_centres[idx]
    wanted = set(obs_eids)

    odb = openOdb(odb_path, readOnly=True)
    try:
        step = odb.steps[odb.steps.keys()[-1]]
        rows = []
        seen = set()
        for name, region in step.historyRegions.items():
            if not name.startswith('Element '):
                continue
            try:
                eid = int(name.split('.')[-1].split()[0])
            except (ValueError, IndexError):
                continue
            if eid not in wanted or 'E33' not in region.historyOutputs:
                continue
            seen.add(eid)
            z = z_of_eid[eid]
            for (t, e33) in region.historyOutputs['E33'].data:
                rows.append((t, eid, z, e33))
    finally:
        odb.close()
    return rows, len(wanted - seen)


def write_csv(path, header, rows):
    # 'wb' on Python 2, 'w' + newline='' on Python 3 — csv needs both.
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
        print('Usage: abaqus python extract_strain.py <job>.odb')
        return 1
    odb_path = sys.argv[1]
    if not os.path.isfile(odb_path):
        print('ERROR: no such ODB: ' + odb_path)
        return 1

    base, out_dir, meta_path = find_meta(odb_path)
    meta = json.load(open(meta_path))

    z_nodes = meta['Z_NODES']
    nz = meta['NZ']
    z_centres = [0.5 * (z_nodes[k] + z_nodes[k + 1]) for k in range(nz)]
    obs = meta['observation']
    obs_eids = obs['OBS_ELEM_IDS']

    print('case %s: observation column at x = %.4f m, %d elements along z'
          % (meta.get('case_id', base), obs['observation_x'], len(obs_eids)))

    rows, missing = collect(odb_path, obs_eids, z_centres)
    if missing:
        print('WARNING: %d observation elements had no E33 history output. '
              'Check that the deck requests E33 on ELSET=OBS_ELEMS.' % missing)
    if not rows:
        print('ERROR: no E33 history collected — nothing to write.')
        return 2
    print('collected %d samples' % len(rows))

    write_csv(os.path.join(out_dir, base + '_e33_history.csv'),
              ['time', 'element_id', 'z_centre', 'E33'], rows)

    centre_eid = obs_eids[len(z_centres) // 2]
    centre = [(t, e33) for (t, eid, z, e33) in rows if eid == centre_eid]
    centre.sort()
    write_csv(os.path.join(out_dir, base + '_e33_centre.csv'),
              ['time', 'E33'], centre)
    return 0


if __name__ == '__main__':
    sys.exit(main())

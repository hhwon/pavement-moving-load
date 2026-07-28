# CLAUDE.md

The agent guide for this repository is **[AGENTS.md](AGENTS.md)** — read that first.

Quick reminders:

- **This repository is public.** No personal paths, hostnames, IPs, usernames or
  e-mail addresses. `tests/test_build.py` scans every file for them, itself
  included; `example.org` is the accepted placeholder.
- `export PYTHONPATH=$PWD:$PYTHONPATH`, then
  `python -m movingload.submit --structure semi_rigid --motion vertical --build-only`.
  Deck generation needs neither Abaqus nor a network.
- Run `python tests/test_build.py` (32 checks, seconds) after any change.
- `movingload/extract_strain.py` runs under Abaqus' Python 2.7 — no f-strings,
  no `pathlib`, no annotations.
- `movingload/cluster.py` has **no default SSH target** and must keep raising
  `ClusterNotConfigured` instead of guessing.
- Materials are isotropic linear elastic throughout; there is no bimodular
  behaviour in this repository.

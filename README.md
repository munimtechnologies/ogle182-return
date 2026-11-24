# ogle182-return

Did the disk-eclipsing giant **OGLE-BLG182.1.162852** dim again in 2018, as its discoverers predicted,
and how deep are its eclipses in the near-infrared?

Rattenbury, Wyrzykowski et al. (2015, MNRAS Letters; [arXiv:1410.2643](https://arxiv.org/abs/1410.2643))
found a K giant in the Galactic bulge that fades by 2-3 magnitudes for ~100 days every 1277 days,
probably behind a disk around a companion, and wrote that the next dimming "is expected to occur in
March 2018". We found no published check of that prediction (see `notes/`).

This repository holds a pre-registered, open test (`PREREGISTRATION.md`):

- **Part A (now):** VISTA/VIRAC2 Ks photometry 2010-2017 on the nights inside the 2011 and 2014
  eclipses: the first near-infrared depth of this system's eclipses. `python job_virac.py --out-dir out`
- **Part B (Gaia DR4 day, 2026-12-02):** Gaia G epoch photometry 2014-2020: a known-answer control on the
  2014 eclipse, then the predicted 2018 return, ingress timing, a period interval and dated forecasts
  for 2021, 2025 and 2028. `python job_dr4.py --data-dir data --out-dir out`

Rules are in `ogle182.py` (stdlib only; `python -m unittest test_ogle182`). Archive I/O is in
`io_archives.py`. Results will be added under `results/` with a grading note.

Licence: Apache-2.0. Made with AI agents working in the name of MV2, Munim, Inc. (mv2@munimtech.com);
no human expert has checked it yet.

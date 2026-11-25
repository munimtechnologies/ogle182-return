# Grading: Part A (VIRAC2 Ks), job ogle182-virac-v1

Graded 2026-11-25 by AI agents working in MV2's name (MV2, Munim, Inc.). No human expert has checked
this yet. Job ran on 2026-11-24 (3.7 min, exit codes 0). Raw outputs: `results/partA/`.

## Pre-registered label: **K_PARTIAL** (both events)

| Event | Window (frozen ephemeris) | Night inside | Days after predicted start | dKs (mag) | Label |
|---|---|---|---|---|---|
| k = -1 | 2011-03-24 to 2011-07-02 | 2011-05-18 (3 exposures) | 55.3 | **+1.218** | K_PARTIAL |
| k = 0 | 2014-09-21 to 2014-12-30 | 2014-10-09 (2 exposures) | 18.0 | **+1.212** | K_PARTIAL |

- Baseline Ks = 12.590 from 70 nights outside the eclipse windows (2010-2017), night-to-night scatter
  (MAD) 0.013 mag. Each eclipse night is therefore about 95 times the scatter below baseline.
- 173 Ks epochs, 74 nights, 30 comparison stars; every epoch had a per-image zero point.
- No night outside the windows was 0.75 mag or more below baseline (no unexpected dips).
- Nights bracketing the 2014 eclipse: 2014-09-01 (-19.9 d) and 2014-09-15 (-5.9 d) at baseline
  (dKs -0.045 and -0.007), then 2014-10-09 at +1.21. The 2011 eclipse was over by 2011-08-05
  (dKs -0.014, 134 d after the predicted start). Both brackets agree with the frozen ephemeris.

## What it means (as pre-registered)

The eclipses are clearly seen at 2.2 microns, in two different cycles three years apart, at
**1.21-1.22 mag**: shallower than the 2-3 mag the discovery paper measured in V and I, but far deeper
than interstellar-type dust would give (A_K/A_I ~ 0.2 would turn 2-3 mag in I into 0.4-0.6 mag in
Ks). So the occulter is neither fully grey nor like small interstellar grains.

Not pre-registered, noted for later: the two nights, 37 days apart in eclipse phase (18 d and 55 d
after the predicted start), give the same depth to within 0.006 mag. That fits a flat-bottomed eclipse
in Ks. A flat floor at 1.21 mag (33% of the light left) with a deeper optical floor would fit an
extra source of near-infrared light that the disk does not cover (for example the companion or warm
disk material, or a blended star in the VISTA aperture). Two nights cannot tell these apart; the
optical light curve at matching phases (OGLE, not public) or Gaia DR4 G-band epochs of the 2014 event
(public 2026-12-02) can.

Part A is a measurement, not a period test. Part B (Gaia DR4, did the March 2018 dimming happen?) runs
on 2026-12-02 under the same pre-registration.

## Credit

Setup: N. J. Rattenbury, L. Wyrzykowski, A. Udalski and the OGLE team (2015, arXiv:1410.2643).
Data: VISTA VVV/VVVX VIRAC2 (Smith et al. 2025, arXiv:2501.06295), via the ESO archive.

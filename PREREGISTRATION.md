# Pre-registration: did OGLE-BLG182.1.162852 dim again in 2018, and how deep are its eclipses in the infrared?

Frozen 2026-11-24 (UTC), before any brightness measurement of the target after the discovery paper was
retrieved by us. Code: `ogle182.py` (rules), `job_virac.py` (Part A), `job_dr4.py` (Part B), tests in
`test_ogle182.py`. Licence Apache-2.0. Made with AI agents working in MV2's name (MV2, Munim, Inc.);
no human expert has checked this design yet.

## 1. The near miss

Rattenbury, Wyrzykowski, Kostrzewa-Rutkowska, Udalski, Kozlowski, Szymanski, Pietrzynski, Soszynski,
Poleski, Ulaczyk, Skowron, Pietrukowicz, Mroz and Skowron (2015, MNRAS Letters; arXiv:1410.2643) found
that the bulge K giant OGLE-BLG182.1.162852 (17:59:03.54 -30:49:06.0) dims by 2-3 mag in I for about
100 days every 1277 days, most likely because a disk around an unseen companion passes in front of it.
In late 2014 a dimming was under way, "predicted to end on December 30th, 2014", and they wrote:
"The next dimming episode for OGLE-BLG182.1.162852 is expected to occur in March 2018." They asked for
multi-band photometry. They also saw that the depth changes from event to event (3 mag, then 2 mag)
and suggested the disk might precess.

As far as our searches show (notes/scoop-check-2026-11-24.md), no paper has reported whether the 2018
dimming happened, nor measured the eclipse depth in the near-infrared. The one 2026 paper that lists the
system (Hu, Zhu, Wang & Wang, arXiv:2601.16828, Table 2) quotes the 1277-day period from the discovery
paper and does not analyse new data for it. The OGLE-IV light curve is not public.

Two public data sets reach the eclipses:
- **VIRAC2** (Smith et al. 2025, arXiv:2501.06295; ESO table safcat.VVVX_VIRAC_V2_LC): VISTA Ks photometry
  of the target (VIRAC2 source 13710333009739, 0.04" from the OGLE position) from 2010-07 to 2017-03.
- **Gaia DR4** (public 2026-12-02): G-band epoch photometry from 2014-07-25 to 2020-01, which covers
  the 2014 dimming (known answer) and the predicted 2018 dimming. Gaia DR3 source 4044166116946029568
  (0.05" away; DR3 has no epoch photometry for it).

## 2. Frozen ephemeris

P = 1277 d and duration 100 d (Table 1 of the paper); anchor: the 2014 event ends on 2014-12-30
(MJD 57021), so it starts at MJD 56921 (2014-09-21). Event k runs from 56921 + 1277k for 100 days:
k=-1 2011-03-24, k=0 2014-09-21, k=1 2018-03-21, k=2 2021-09, k=3 2025-03-18, k=4 2028-09.
Each event has a **wide** window (start - 20 d .. end + 20 d) and a **core** window
(start + 20 d .. end - 20 d). The 20-day tolerance covers the unstated period and anchor uncertainty.

## 3. Part A: VIRAC2 Ks (runs now)

- Data: VIRAC2 sources within a 60" box; Ks light curves of the target and up to 30 comparison stars
  (>= 5" away, 11.5 <= Ks <= 14.0, >= 120 detections, not duplicate, Ks MAD <= 0.05, smallest MAD first).
- Cleaning: Ks rows with finite mag, emag < 0.2, ambiguous_match = 0; phot_flag = 0 also required if that
  keeps >= 50 % of rows (else ignored, logged).
- Zero point per image (catid): median over comparison stars of (mag - star median), needs >= 5 stars;
  otherwise no correction (counted).
- Nights: epochs within 0.5 d merged (median). Baseline: median of nights outside the wide windows of
  k=-1 and k=0 (needs >= 10). dKs = night - baseline (positive = fainter).
- "Inside an eclipse": a night 10-90 days after the predicted start of k=-1 or k=0.
- Labels (deepest inside night): **K_DEEP** dKs >= 1.5; **K_PARTIAL** 0.75-1.5; **K_WEAK** 0.30-0.75;
  **K_FLAT** <= 0.30; **K_NO_COVERAGE** no inside night. Per event and overall (deepest).
- Meaning: the paper found the eclipse grey from V to I (dV ~ dI, 2-3 mag). K_DEEP means it stays near
  grey out to 2.2 microns (large grains or an optically thick screen); K_WEAK or K_FLAT on a night well
  inside a predicted eclipse means either strongly chromatic dimming (small grains; interstellar-type
  dust has A_K/A_I ~ 0.2) or that the eclipse was not under way, i.e. the ephemeris is off. Part A cannot
  tell those two apart alone; we will say so.
- Disclosed look: before freezing we listed the epoch TIMES of the target's VIRAC2 light curve (no
  magnitudes). Ks epochs inside the windows: 2011-05-18 (55 d after the predicted start, 3 exposures)
  and 2014-10-09 (18 d after, 2 exposures); 2014-09-01 and 2014-09-15 lie just before the 2014 start.
  So Part A rests on two nights; it is a measurement, not a period test.
- Also reported: nights with dKs >= 0.75 outside all windows (unexpected dips).

## 4. Part B: Gaia DR4 G epoch photometry (DR4 day, 2026-12-02 or later)

- Data: DataLink epoch photometry of the target (DR4 id by a 1" cone search, else the DR3 id) and of a
  scan tracer, the G = 16.5 neighbour Gaia DR3 4044166121093245952 (2.9" away), whose transits mark when
  Gaia scanned the field. G-band transits only; photometry-rejected and non-positive fluxes dropped;
  variability-rejected transits kept (eclipse points can be flagged as outliers).
- Baseline: median flux of transits outside the wide windows of k=0 and k=1 (needs >= 10).
  dG = -2.5 log10(flux / baseline). DIM: dG >= 0.75; BRIGHT: dG <= 0.30. Visits: transits within 2 d.
- A tracer visit with no target transit within 1 day counts as a **missed visit** (in eclipse the
  target, G ~ 17.2 out of eclipse, may fall below Gaia's detection limit in this crowded field).
- **Control (2014 event, known answer):** PASS if >= 2 DIM transits in the k=0 wide window; FAIL if
  none and >= 2 BRIGHT visits in the core; else NO_COVERAGE. If FAIL, the main label is reported but
  not interpreted (pipeline or identification suspect).
- **Main (2018 prediction):**
  - RETURNED_ON_SCHEDULE: >= 2 DIM transits in the k=1 wide window (2018-03-01 .. 2018-07-19);
  - else RETURNED_SHIFTED: an episode of >= 2 DIM transits (linked within 100 d) outside the windows
    starting within 200 d of the predicted window (shift reported);
  - else NOT_RETURNED: no DIM transit in the wide window, >= 2 BRIGHT visits in the core
    (2018-04-10 .. 2018-06-09) and no missed core visit;
  - else DROPOUT_ONLY: no DIM transit, >= 2 missed core visits and no BRIGHT core visit (consistent with
    fading below detection; not counted as a confirmation);
  - else INCONCLUSIVE.
- Ingress brackets (last BRIGHT transit before the first DIM transit, within 60 d) for 2014 and 2018
  give a period interval from Gaia alone (EMPTY if they contradict a strict period); with it, dated
  ingress windows for 2021, 2025 and 2028 are written as a forecast anyone can check.
- Off-window DIM episodes are listed whatever the main label.

## 5. What would count as news

RETURNED_ON_SCHEDULE: the first confirmation, from public data, of the 2018 return that the discovery
paper predicted, with a measured depth and a period interval good enough to time the 2028 eclipse.
NOT_RETURNED or a large shift: the system is not strictly periodic at 1277 d, which challenges the
eclipsing-binary picture. K_DEEP: the occulter is grey out to 2.2 microns.

## 6. Deviations

None yet. Any change after this commit will be logged here with its date and reason.

# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""OGLE-BLG182.1.162852 return test: frozen rules (PREREGISTRATION.md, 2026-11-24). Stdlib only.

Question: did the disk-eclipsing giant OGLE-BLG182.1.162852 dim again in 2018, as Rattenbury,
Wyrzykowski et al. (2015, MNRAS 447, L31; arXiv:1410.2643) predicted ("The next dimming episode ...
is expected to occur in March 2018"), and how deep are its eclipses in the near-infrared?

Part A (VIRAC2, VISTA Ks, 2010-2017): Ks depth on the nights that fall inside the 2011 and 2014
  eclipses that OGLE saw (the first near-infrared depth of this system's eclipses).
Part B (Gaia DR4 epoch photometry, 2014-07 .. 2020-01, public 2026-12-02): known-answer control on
  the 2014 eclipse, then the predicted 2018 return; ingress brackets, period interval, forecasts.
Made with AI agents working in MV2's name (MV2, Munim, Inc.).
"""
import math
from datetime import date, timedelta

VERSION = "0.1.0"

# Target (Rattenbury et al. 2015): 17:59:03.54 -30:49:06.0 (J2000)
RA, DEC = 269.764750, -30.818333
VIRAC_ID = 13710333009739            # VIRAC2 source at 269.764762 -30.818336 (0.04" away)
GAIA_DR3_ID = 4044166116946029568    # Gaia DR3 source at 269.764764 -30.818339 (0.05" away)
GAIA_TRACER_ID = 4044166121093245952  # G=16.5 neighbour 2.9" away: its transits mark Gaia's scans

_MJD0 = date(1858, 11, 17)


def mjd(y, m, d):
    return (date(y, m, d) - _MJD0).days


def date_str(t):
    return (_MJD0 + timedelta(days=int(math.floor(t)))).isoformat()


# --- the published ephemeris (frozen) ---------------------------------------------------------
P_DAYS = 1277.0                 # "separated by 1277 days"
DURATION = 100.0                # "~100 days" (Table 1)
END0 = mjd(2014, 12, 30)        # the 2014 dimming "predicted to end on December 30th, 2014" (57021)
START0 = END0 - DURATION        # 56921 = 2014-09-21
TOL = 20.0                      # timing tolerance on each side (period and anchor uncertainty)
EVENT_YEARS = {-1: 2011, 0: 2014, 1: 2018, 2: 2021, 3: 2025, 4: 2028}

# --- frozen thresholds -------------------------------------------------------------------------
DIM_MAG = 0.75        # >= 0.75 mag fainter than baseline (published depth 2-3 mag in I)
BRIGHT_MAG = 0.30     # within 0.30 mag of baseline (out-of-eclipse spots: ~0.06 mag)
DEEP_K = 1.5          # Ks depth that counts as near-grey against the published 2-3 mag in I
PHASE_IN = (10.0, 90.0)   # days after predicted start that count as "inside the eclipse"
NIGHT_GAP = 0.5       # VIRAC epochs closer than this (days) form one night
VISIT_GAP = 2.0       # Gaia transits closer than this (days) form one visit
MIN_BASE = 10         # baseline needs >= 10 nights / transits outside all event windows
MIN_COMPS = 5         # per-epoch zero point needs >= 5 comparison stars
MAX_BRACKET = 60.0    # an ingress bracket is used only if it is <= 60 days wide
SHIFT_SEARCH = 200.0  # a DIM episode within +-200 d of the 2018 window counts as a shifted return
EPISODE_LINK = 100.0  # DIM transits closer than this (days) belong to one episode


def window(k):
    s = START0 + k * P_DAYS
    e = s + DURATION
    return {"k": k, "year": EVENT_YEARS.get(k), "start": s, "end": e,
            "date_start": date_str(s), "date_end": date_str(e),
            "wide": (s - TOL, e + TOL), "core": (s + TOL, e - TOL)}


def median(xs):
    s = sorted(xs)
    n = len(s)
    if n == 0:
        return float("nan")
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def in_range(t, ab):
    return ab[0] <= t <= ab[1]


def in_any_wide(t, ks):
    return any(in_range(t, window(k)["wide"]) for k in ks)


def group(points, gap):
    """points: list of (t, value) -> groups of consecutive points with spacing <= gap.
    Each group: t (mean), t_lo, t_hi, n, value (median)."""
    pts = sorted(points)
    out = []
    cur = []
    for p in pts:
        if cur and p[0] - cur[-1][0] > gap:
            out.append(cur)
            cur = []
        cur.append(p)
    if cur:
        out.append(cur)
    res = []
    for g in out:
        ts = [p[0] for p in g]
        res.append({"t": sum(ts) / len(ts), "t_lo": ts[0], "t_hi": ts[-1], "n": len(g),
                    "value": median([p[1] for p in g])})
    return res


def state(dm):
    if dm >= DIM_MAG:
        return "DIM"
    if dm <= BRIGHT_MAG:
        return "BRIGHT"
    return "MID"


# ============================ Part A: VIRAC2 Ks ===============================================
def select_comps(sources, target_id=VIRAC_ID, ra0=RA, dec0=DEC, max_n=30):
    """Comparison stars: not the target, >= 5" from it, 11.5 <= Ks <= 14.0, >= 120 Ks detections,
    not a duplicate, Ks MAD <= 0.05 mag; the max_n with the smallest MAD."""
    cosd = math.cos(math.radians(dec0))
    keep = []
    for s in sources:
        try:
            sid = int(s["sourceid"])
            ra, de = float(s["ra"]), float(s["de"])
            km = float(s["phot_ks_mean_mag"])
            nd = int(float(s["ks_n_det"]))
            mad = float(s["ks_mad"])
            dup = int(float(s.get("duplicate") or 0))
        except (KeyError, TypeError, ValueError):
            continue
        sep = 3600.0 * math.hypot((ra - ra0) * cosd, de - dec0)
        if sid == target_id or sep < 5.0 or dup != 0:
            continue
        if 11.5 <= km <= 14.0 and nd >= 120 and mad <= 0.05:
            keep.append((mad, sid))
    keep.sort()
    return [sid for _, sid in keep[:max_n]]


def clean_lc(rows):
    """rows: dicts with sourceid, catid, mjdobs, filter, mag, emag, phot_flag, ambiguous_match.
    Keeps Ks rows with a finite mag, emag < 0.2 and ambiguous_match == 0. phot_flag == 0 is also
    required if that keeps >= 50% of the rows (logged in the result), else phot_flag is ignored."""
    base = []
    for r in rows:
        try:
            if (r.get("filter") or "").strip() != "Ks":
                continue
            mag = float(r["mag"])
            emag = float(r["emag"])
            if not (mag == mag) or not (emag == emag) or emag >= 0.2:
                continue
            if int(float(r.get("ambiguous_match") or 0)) != 0:
                continue
            flag = int(float(r.get("phot_flag") or 0))
            base.append({"sourceid": int(r["sourceid"]), "catid": int(r["catid"]),
                         "t": float(r["mjdobs"]), "mag": mag, "emag": emag, "flag": flag})
        except (KeyError, TypeError, ValueError):
            continue
    good = [b for b in base if b["flag"] == 0]
    used = len(good) >= 0.5 * len(base) and len(good) > 0
    return (good if used else base), {"n_rows": len(base), "n_flag0": len(good), "phot_flag_used": used}


def zero_points(comp_rows, min_comps=MIN_COMPS):
    """Per-image (catid) offset = median over comparison stars of (mag - that star's median)."""
    by_star = {}
    for r in comp_rows:
        by_star.setdefault(r["sourceid"], []).append(r)
    resid = {}
    for sid, rs in by_star.items():
        m = median([r["mag"] for r in rs])
        for r in rs:
            resid.setdefault(r["catid"], []).append(r["mag"] - m)
    return {c: median(v) for c, v in resid.items() if len(v) >= min_comps}


def virac_grade(target_rows, comp_rows, ks=(-1, 0)):
    zp = zero_points(comp_rows)
    corr = []
    n_nozp = 0
    for r in target_rows:
        off = zp.get(r["catid"])
        if off is None:
            n_nozp += 1
            off = 0.0
        corr.append((r["t"], r["mag"] - off))
    nights = group(corr, NIGHT_GAP)
    base_n = [n["value"] for n in nights if not in_any_wide(n["t"], ks)]
    res = {"n_epochs": len(corr), "n_epochs_without_zp": n_nozp, "n_nights": len(nights),
           "n_baseline_nights": len(base_n)}
    if len(base_n) < MIN_BASE:
        res["label"] = "K_NO_BASELINE"
        return res, []
    base = median(base_n)
    res["baseline_ks"] = base
    out_nights = []
    for n in nights:
        dm = n["value"] - base
        ev = None
        phase = None
        for k in ks:
            w = window(k)
            if in_range(n["t"], w["wide"]):
                ev, phase = k, n["t"] - w["start"]
        out_nights.append({"t": n["t"], "date": date_str(n["t"]), "n": n["n"], "ks": n["value"],
                           "dm": dm, "state": state(dm), "event": ev, "phase": phase})
    events = {}
    inside_all = []
    for k in ks:
        w = window(k)
        ins = [x for x in out_nights if x["event"] == k and x["phase"] is not None
               and PHASE_IN[0] <= x["phase"] <= PHASE_IN[1]]
        inside_all += ins
        events[str(k)] = {"year": w["year"], "date_start": w["date_start"], "date_end": w["date_end"],
                          "n_nights_wide": sum(1 for x in out_nights if x["event"] == k),
                          "n_nights_inside": len(ins),
                          "nights_inside": [{"date": x["date"], "phase_d": round(x["phase"], 1),
                                             "dm": round(x["dm"], 3)} for x in ins],
                          "max_dm_inside": max((x["dm"] for x in ins), default=None),
                          "label": k_label(max((x["dm"] for x in ins), default=None))}
    res["events"] = events
    res["max_dm_inside"] = max((x["dm"] for x in inside_all), default=None)
    res["label"] = k_label(res["max_dm_inside"])
    res["unexpected_dim_nights"] = [{"date": x["date"], "dm": round(x["dm"], 3)} for x in out_nights
                                    if x["event"] is None and x["state"] == "DIM"]
    res["baseline_scatter_mad"] = median([abs(x["dm"]) for x in out_nights if x["event"] is None])
    return res, out_nights


def k_label(dm):
    if dm is None:
        return "K_NO_COVERAGE"
    if dm >= DEEP_K:
        return "K_DEEP"
    if dm >= DIM_MAG:
        return "K_PARTIAL"
    if dm > BRIGHT_MAG:
        return "K_WEAK"
    return "K_FLAT"


# ============================ Part B: Gaia DR4 G epoch photometry =============================
GAIA_FLUX_ZP = 25.6874


def gaia_time_to_mjd(v):
    v = float(v)
    if v > 2400000.0:
        return v - 2400000.5
    if v > 40000.0:
        return v
    return v + 55197.0   # Gaia time = BJD(TCB) - 2455197.5


def _truthy(v):
    if v is None:
        return False
    if isinstance(v, bool):
        return v
    s = str(v).strip().lower()
    return s in ("true", "t", "1", "yes")


def _num(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if x == x and abs(x) != float("inf") else None


def _arr(v):
    if v is None:
        return []
    if isinstance(v, (list, tuple)):
        return list(v)
    s = str(v).strip().strip("()[]")
    if not s:
        return []
    return [p.strip() for p in s.replace(";", ",").split(",")]


def parse_gaia_rows(rows, source_id=None):
    """Gaia epoch photometry rows -> sorted [(mjd, flux)] for the G band.
    Accepts (a) one row per transit with columns band/time/flux|mag(/rejected_by_photometry), and
    (b) one row per source with arrays g_transit_time and g_transit_flux|g_transit_mag.
    Rows of other sources (if a source_id column exists) are skipped. Photometry-rejected or
    non-positive fluxes are dropped; variability-rejected points are KEPT (eclipse points)."""
    out = []
    for r in rows:
        if source_id is not None and r.get("source_id") not in (None, ""):
            try:
                if int(float(r["source_id"])) != int(source_id):
                    continue
            except (TypeError, ValueError):
                pass
        if "g_transit_time" in r:
            ts = _arr(r.get("g_transit_time"))
            fs = _arr(r.get("g_transit_flux")) if r.get("g_transit_flux") is not None else []
            ms = _arr(r.get("g_transit_mag")) if r.get("g_transit_mag") is not None else []
            rej = _arr(r.get("rejected_by_photometry")) if r.get("rejected_by_photometry") is not None else []
            for i, t in enumerate(ts):
                tv = _num(t)
                if tv is None:
                    continue
                f = _num(fs[i]) if i < len(fs) else None
                if f is None and i < len(ms) and _num(ms[i]) is not None:
                    f = 10 ** (-0.4 * (_num(ms[i]) - GAIA_FLUX_ZP))
                if f is None or f <= 0:
                    continue
                if i < len(rej) and _truthy(rej[i]):
                    continue
                out.append((gaia_time_to_mjd(tv), f))
            continue
        band = (r.get("band") or "").strip().upper()
        if band != "G":
            continue
        tv = _num(r.get("time"))
        if tv is None:
            continue
        if _truthy(r.get("rejected_by_photometry")):
            continue
        f = _num(r.get("flux"))
        if f is None and _num(r.get("mag")) is not None:
            f = 10 ** (-0.4 * (_num(r.get("mag")) - GAIA_FLUX_ZP))
        if f is None or f <= 0:
            continue
        out.append((gaia_time_to_mjd(tv), f))
    out.sort()
    return out


def to_dm(transits, base):
    return [(t, -2.5 * math.log10(f / base)) for t, f in transits]


def missed_visits(target, tracer, gap=VISIT_GAP):
    """Tracer visits (Gaia scans over the field) in which the target has no transit within 1 day."""
    tv = group([(t, 0.0) for t, _ in tracer], gap)
    tt = [t for t, _ in target]
    out = []
    for v in tv:
        hit = any(v["t_lo"] - 1.0 <= t <= v["t_hi"] + 1.0 for t in tt)
        out.append({"t": v["t"], "date": date_str(v["t"]), "detected": hit})
    return out


def event_counts(pts, visits, tracer_visits, k):
    w = window(k)
    n_dim = sum(1 for t, dm in pts if in_range(t, w["wide"]) and dm >= DIM_MAG)
    n_wide = sum(1 for t, dm in pts if in_range(t, w["wide"]))
    core_v = [v for v in visits if in_range(v["t"], w["core"])]
    n_bright_core = sum(1 for v in core_v if v["value"] <= BRIGHT_MAG)
    n_missed_core = sum(1 for v in tracer_visits if in_range(v["t"], w["core"]) and not v["detected"])
    return {"year": w["year"], "date_start": w["date_start"], "date_end": w["date_end"],
            "n_transits_wide": n_wide, "n_dim_transits": n_dim, "n_core_visits": len(core_v),
            "n_bright_core_visits": n_bright_core, "n_missed_core_visits": n_missed_core,
            "max_dm_wide": max((dm for t, dm in pts if in_range(t, w["wide"])), default=None)}


def dim_episodes(pts, ks=(0, 1)):
    """Episodes of >= 2 DIM transits (linked within EPISODE_LINK days) outside all wide windows."""
    dims = sorted(t for t, dm in pts if dm >= DIM_MAG and not in_any_wide(t, ks))
    eps = []
    cur = []
    for t in dims:
        if cur and t - cur[-1] > EPISODE_LINK:
            eps.append(cur)
            cur = []
        cur.append(t)
    if cur:
        eps.append(cur)
    return [{"t_first": e[0], "t_last": e[-1], "date_first": date_str(e[0]), "date_last": date_str(e[-1]),
             "n_dim": len(e)} for e in eps if len(e) >= 2]


def ingress_bracket(pts, k):
    """(last BRIGHT transit before the first DIM transit, first DIM transit) for event k, searched
    in [wide_lo - MAX_BRACKET, wide_hi]; None unless both exist and the gap is <= MAX_BRACKET."""
    w = window(k)
    lo, hi = w["wide"][0] - MAX_BRACKET, w["wide"][1]
    sel = [(t, dm) for t, dm in pts if lo <= t <= hi]
    first_dim = next((t for t, dm in sel if dm >= DIM_MAG), None)
    if first_dim is None:
        return None
    before = [t for t, dm in sel if t < first_dim and dm <= BRIGHT_MAG]
    if not before:
        return None
    b = max(before)
    if first_dim - b > MAX_BRACKET:
        return None
    return (b, first_dim)


def period_interval(brackets):
    """brackets: {k: (lo, hi)} ingress brackets -> (P_lo, P_hi) consistent with all pairs, or None
    if fewer than two, or 'EMPTY' if they contradict a strict period."""
    ks = sorted(brackets)
    if len(ks) < 2:
        return None
    plo, phi = -float("inf"), float("inf")
    for i in range(len(ks)):
        for j in range(i + 1, len(ks)):
            a, b = brackets[ks[i]], brackets[ks[j]]
            n = ks[j] - ks[i]
            plo = max(plo, (b[0] - a[1]) / n)
            phi = min(phi, (b[1] - a[0]) / n)
    if plo > phi:
        return "EMPTY"
    return (plo, phi)


def forecast(brackets, pint, ks=(2, 3, 4)):
    """Ingress windows for later events from the latest bracket and the period interval."""
    if not isinstance(pint, tuple) or not brackets:
        return []
    k0 = max(brackets)
    lo, hi = brackets[k0]
    out = []
    for k in ks:
        n = k - k0
        a, b = lo + n * pint[0], hi + n * pint[1]
        out.append({"k": k, "year": EVENT_YEARS.get(k), "ingress_from": date_str(a), "ingress_to": date_str(b),
                    "ingress_from_mjd": a, "ingress_to_mjd": b})
    return out


def gaia_grade(target, tracer):
    """target, tracer: sorted [(mjd, flux)] G transits. Returns the frozen labels and details."""
    ks = (0, 1)
    base_f = [f for t, f in target if not in_any_wide(t, ks)]
    res = {"n_transits": len(target), "n_baseline_transits": len(base_f), "n_tracer_transits": len(tracer)}
    if len(base_f) < MIN_BASE:
        res.update({"control": "NO_BASELINE", "main": "NO_BASELINE"})
        return res, []
    base = median(base_f)
    res["baseline_flux"] = base
    pts = to_dm(target, base)
    visits = group(pts, VISIT_GAP)
    tvis = missed_visits(target, tracer) if tracer else []
    c0 = event_counts(pts, visits, tvis, 0)
    c1 = event_counts(pts, visits, tvis, 1)
    eps = dim_episodes(pts, ks)
    res.update({"event_2014": c0, "event_2018": c1, "off_window_dim_episodes": eps,
                "n_tracer_visits": len(tvis), "n_tracer_visits_missed": sum(1 for v in tvis if not v["detected"])})
    # control: the 2014 eclipse that OGLE was watching when the paper was written
    if c0["n_dim_transits"] >= 2:
        res["control"] = "PASS"
    elif c0["n_dim_transits"] == 0 and c0["n_bright_core_visits"] >= 2:
        res["control"] = "FAIL"
    else:
        res["control"] = "NO_COVERAGE"
    # main: the predicted 2018 return
    w1 = window(1)
    shifted = [e for e in eps if w1["start"] - SHIFT_SEARCH <= e["t_first"] <= w1["end"] + SHIFT_SEARCH]
    if c1["n_dim_transits"] >= 2:
        main = "RETURNED_ON_SCHEDULE"
    elif shifted:
        main = "RETURNED_SHIFTED"
        res["shift_days"] = shifted[0]["t_first"] - w1["start"]
    elif c1["n_dim_transits"] == 0 and c1["n_bright_core_visits"] >= 2 and c1["n_missed_core_visits"] == 0:
        main = "NOT_RETURNED"
    elif c1["n_dim_transits"] == 0 and c1["n_missed_core_visits"] >= 2 and c1["n_bright_core_visits"] == 0:
        main = "DROPOUT_ONLY"
    else:
        main = "INCONCLUSIVE"
    res["main"] = main
    br = {}
    for k in ks:
        b = ingress_bracket(pts, k)
        if b:
            br[k] = b
    res["ingress_brackets"] = {str(k): {"from": date_str(v[0]), "to": date_str(v[1]), "from_mjd": v[0],
                                       "to_mjd": v[1]} for k, v in br.items()}
    pint = period_interval(br)
    res["period_interval_days"] = list(pint) if isinstance(pint, tuple) else pint
    res["forecast"] = forecast(br, pint)
    return res, pts

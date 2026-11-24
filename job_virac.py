# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""Part A: VIRAC2 Ks depth of OGLE-BLG182.1.162852's 2011 and 2014 eclipses (frozen rules).

python job_virac.py --out-dir /out
Made with AI agents working in MV2's name (MV2, Munim, Inc.).
"""
import argparse
import csv
import json
import math
import os
import time

import io_archives as io_a
import ogle182 as o


def log(msg):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), msg), flush=True)


def write_csv(path, rows, cols):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in rows:
            w.writerow([r.get(c) for c in cols])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    cosd = math.cos(math.radians(o.DEC))
    srcs = io_a.virac_sources_box(o.RA, o.DEC, 60.0 / 3600.0, cosd, log=log)
    log("VIRAC2 sources in the 60-arcsec box: %d" % len(srcs))
    write_csv(os.path.join(a.out_dir, "virac_sources.csv"), srcs,
              ["sourceid", "ra", "de", "phot_ks_mean_mag", "ks_n_det", "ks_mad", "duplicate"])
    comps = o.select_comps(srcs)
    log("comparison stars: %d" % len(comps))
    rows = io_a.virac_lightcurves([o.VIRAC_ID] + comps, log=log)
    log("light-curve rows: %d" % len(rows))
    write_csv(os.path.join(a.out_dir, "virac_lc_raw.csv"), rows,
              ["sourceid", "catid", "mjdobs", "filter", "mag", "emag", "phot_flag", "ambiguous_match", "chi",
               "seeing"])
    tgt_raw = [r for r in rows if str(r.get("sourceid")).strip() == str(o.VIRAC_ID)]
    cmp_raw = [r for r in rows if str(r.get("sourceid")).strip() != str(o.VIRAC_ID)]
    tgt, tmeta = o.clean_lc(tgt_raw)
    cmp_, cmeta = o.clean_lc(cmp_raw)
    res, nights = o.virac_grade(tgt, cmp_)
    res.update({"version": o.VERSION, "part": "A (VIRAC2 Ks)", "virac_id": o.VIRAC_ID,
                "n_comps": len(comps), "comps": comps, "target_clean": tmeta, "comps_clean": cmeta,
                "windows": {str(k): o.window(k) for k in (-1, 0, 1, 2, 3, 4)},
                "thresholds": {"DIM_MAG": o.DIM_MAG, "BRIGHT_MAG": o.BRIGHT_MAG, "DEEP_K": o.DEEP_K,
                               "PHASE_IN": o.PHASE_IN, "TOL": o.TOL}})
    with open(os.path.join(a.out_dir, "results_virac.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, indent=1, default=str)
    write_csv(os.path.join(a.out_dir, "virac_nights.csv"), nights,
              ["t", "date", "n", "ks", "dm", "state", "event", "phase"])
    chart = [{"mjd": round(n["t"], 3), "year": round(2000 + (n["t"] - 51544.0) / 365.25, 3),
              "dks": round(n["dm"], 3), "in_eclipse_window": n["event"] is not None} for n in nights]
    write_csv(os.path.join(a.out_dir, "chart_k_lightcurve.csv"), chart,
              ["mjd", "year", "dks", "in_eclipse_window"])
    lines = ["# OGLE-BLG182.1.162852, Part A (VIRAC2 Ks)", "",
             "Label: **%s** (max Ks dimming inside the eclipses: %s mag)" % (res.get("label"), res.get("max_dm_inside")),
             "", "Nights: %s; baseline nights: %s; comparison stars: %d" % (res.get("n_nights"),
                                                                          res.get("n_baseline_nights"), len(comps))]
    for k, e in (res.get("events") or {}).items():
        lines.append("- event %s (%s .. %s): %s; nights inside %s" % (e["year"], e["date_start"], e["date_end"],
                                                                      e["label"], e["nights_inside"]))
    lines.append("- unexpected dim nights outside the windows: %s" % res.get("unexpected_dim_nights"))
    with open(os.path.join(a.out_dir, "summary_virac.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    log("VERDICT part A: %s" % res.get("label"))


if __name__ == "__main__":
    main()

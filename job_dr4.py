# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""Part B (run on/after Gaia DR4 day, 2026-12-02): Gaia G epoch photometry of
OGLE-BLG182.1.162852: 2014 control, predicted 2018 return, ingress brackets, period, forecasts.

python job_dr4.py --data-dir /data --out-dir /out [--release "Gaia DR4"] [--source-id N]
Made with AI agents working in MV2's name (MV2, Munim, Inc.).
"""
import argparse
import csv
import json
import os
import time

import io_archives as io_a
import ogle182 as o


def log(msg):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), msg), flush=True)


def resolve_id(release, fallback, log=log):
    """DR4 cone search (1 arcsec, brightest G) if a DR4 source table answers; else the DR3 id."""
    for table in ("gaiadr4.gaia_source", "gaiadr4.gaia_source_lite"):
        try:
            rows = io_a.gaia_tap_csv(
                "SELECT TOP 5 source_id, phot_g_mean_mag FROM %s WHERE 1=CONTAINS(POINT(ra,dec),"
                "CIRCLE(%.6f,%.6f,%.6f)) ORDER BY phot_g_mean_mag" % (table, o.RA, o.DEC, 1.0 / 3600), log=log)
            if rows and rows[0].get("source_id"):
                return int(rows[0]["source_id"]), table
        except Exception as e:
            log("cone search on %s failed: %s" % (table, str(e)[:200]))
    return fallback, "dr3-id-fallback"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--release", default="Gaia DR4")
    ap.add_argument("--source-id", type=int, default=None)
    ap.add_argument("--tracer-id", type=int, default=o.GAIA_TRACER_ID)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    d = os.path.join(a.data_dir, "gaia-ogle182")
    sid, how = (a.source_id, "argument") if a.source_id else resolve_id(a.release, o.GAIA_DR3_ID)
    log("target source_id %s (%s)" % (sid, how))
    rows, url = io_a.gaia_epoch_rows(sid, d, a.release, log=log)
    target = o.parse_gaia_rows(rows, sid)
    trows, turl = io_a.gaia_epoch_rows(a.tracer_id, d, a.release, log=log)
    tracer = o.parse_gaia_rows(trows, a.tracer_id)
    log("G transits: target %d, tracer %d" % (len(target), len(tracer)))
    res, pts = o.gaia_grade(target, tracer)
    res.update({"version": o.VERSION, "part": "B (Gaia epoch photometry)", "release": a.release,
                "source_id": sid, "source_id_from": how, "tracer_id": a.tracer_id, "datalink_url": url,
                "tracer_datalink_url": turl, "windows": {str(k): o.window(k) for k in (0, 1, 2, 3, 4)}})
    if not target:
        res["main"] = res["control"] = "NO_DATA"
    with open(os.path.join(a.out_dir, "results_dr4.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, indent=1, default=str)
    with open(os.path.join(a.out_dir, "gaia_transits.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["mjd", "date", "flux", "dg"])
        dm = dict(pts)
        for t, fl in target:
            w.writerow([round(t, 5), o.date_str(t), fl, round(dm[t], 4) if t in dm else ""])
    with open(os.path.join(a.out_dir, "chart_gaia_lightcurve.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["year", "dg"])
        for t, v in pts:
            w.writerow([round(2000 + (t - 51544.0) / 365.25, 4), round(v, 3)])
    with open(os.path.join(a.out_dir, "summary_dr4.md"), "w", encoding="utf-8") as f:
        f.write("# OGLE-BLG182.1.162852, Part B (%s)\n\ncontrol (2014): **%s**\nmain (2018): **%s**\n\n%s\n"
                % (a.release, res.get("control"), res.get("main"),
                   json.dumps({k: res.get(k) for k in ("event_2014", "event_2018", "ingress_brackets",
                                                         "period_interval_days", "forecast",
                                                         "off_window_dim_episodes")}, indent=1, default=str)))
    log("VERDICT part B: control %s, main %s" % (res.get("control"), res.get("main")))


if __name__ == "__main__":
    main()

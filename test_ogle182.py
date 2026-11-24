# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""Tests for the frozen rules (stdlib only). Made with AI agents working in MV2's name."""
import math
import unittest

import ogle182 as o


class TestEphemeris(unittest.TestCase):
    def test_anchor_dates(self):
        self.assertEqual(o.END0, 57021)
        self.assertEqual(o.START0, 56921)
        self.assertEqual(o.date_str(o.START0), "2014-09-21")
        w1 = o.window(1)
        self.assertEqual(w1["start"], 58198)
        self.assertEqual(w1["date_start"], "2018-03-21")
        self.assertEqual(o.window(-1)["date_start"], "2011-03-24")
        self.assertEqual(o.window(3)["date_start"], "2025-03-18")
        self.assertEqual(w1["wide"], (58178.0, 58318.0))
        self.assertEqual(w1["core"], (58218.0, 58278.0))

    def test_group(self):
        g = o.group([(1.0, 0.1), (1.2, 0.3), (5.0, 1.0)], 0.5)
        self.assertEqual(len(g), 2)
        self.assertEqual(g[0]["n"], 2)
        self.assertAlmostEqual(g[0]["value"], 0.2)
        self.assertAlmostEqual(g[0]["t"], 1.1)

    def test_k_label(self):
        self.assertEqual(o.k_label(None), "K_NO_COVERAGE")
        self.assertEqual(o.k_label(2.0), "K_DEEP")
        self.assertEqual(o.k_label(1.0), "K_PARTIAL")
        self.assertEqual(o.k_label(0.5), "K_WEAK")
        self.assertEqual(o.k_label(0.1), "K_FLAT")

    def test_period_interval(self):
        self.assertIsNone(o.period_interval({0: (0.0, 10.0)}))
        p = o.period_interval({0: (0.0, 10.0), 1: (1270.0, 1280.0)})
        self.assertAlmostEqual(p[0], 1260.0)
        self.assertAlmostEqual(p[1], 1280.0)
        self.assertEqual(o.period_interval({0: (0.0, 10.0), 1: (1270.0, 1280.0), 2: (2600.0, 2610.0)}), "EMPTY")


class TestVirac(unittest.TestCase):
    def test_clean_lc(self):
        rows = [
            {"sourceid": "1", "catid": "10", "mjdobs": "55500.1", "filter": "Ks", "mag": "12.5", "emag": "0.02",
             "phot_flag": "0", "ambiguous_match": "0"},
            {"sourceid": "1", "catid": "11", "mjdobs": "55520.1", "filter": "Ks", "mag": "12.5", "emag": "0.02",
             "phot_flag": "0", "ambiguous_match": "0"},
            {"sourceid": "1", "catid": "12", "mjdobs": "55540.1", "filter": "Ks", "mag": "12.5", "emag": "0.02",
             "phot_flag": "0", "ambiguous_match": "0"},
            {"sourceid": "1", "catid": "13", "mjdobs": "55560.1", "filter": "Ks", "mag": "12.5", "emag": "0.02",
             "phot_flag": "1", "ambiguous_match": "0"},
            {"sourceid": "1", "catid": "14", "mjdobs": "55580.1", "filter": "Ks", "mag": "12.5", "emag": "0.30",
             "phot_flag": "0", "ambiguous_match": "0"},
            {"sourceid": "1", "catid": "15", "mjdobs": "55600.1", "filter": "Ks", "mag": "12.5", "emag": "0.02",
             "phot_flag": "0", "ambiguous_match": "1"},
            {"sourceid": "1", "catid": "16", "mjdobs": "55620.1", "filter": "J", "mag": "13.5", "emag": "0.02",
             "phot_flag": "0", "ambiguous_match": "0"},
        ]
        good, meta = o.clean_lc(rows)
        self.assertEqual(len(good), 3)
        self.assertTrue(meta["phot_flag_used"])
        self.assertEqual(meta["n_rows"], 4)

    def test_select_comps(self):
        cosd = math.cos(math.radians(o.DEC))
        src = [
            {"sourceid": str(o.VIRAC_ID), "ra": str(o.RA), "de": str(o.DEC), "phot_ks_mean_mag": "12.5",
             "ks_n_det": "173", "ks_mad": "0.01", "duplicate": "0"},
            {"sourceid": "2", "ra": str(o.RA), "de": str(o.DEC + 3.0 / 3600), "phot_ks_mean_mag": "13.0",
             "ks_n_det": "170", "ks_mad": "0.01", "duplicate": "0"},
            {"sourceid": "3", "ra": str(o.RA + 0.005), "de": str(o.DEC), "phot_ks_mean_mag": "13.0",
             "ks_n_det": "150", "ks_mad": "0.02", "duplicate": "0"},
            {"sourceid": "4", "ra": str(o.RA - 0.005), "de": str(o.DEC), "phot_ks_mean_mag": "15.0",
             "ks_n_det": "150", "ks_mad": "0.02", "duplicate": "0"},
        ]
        self.assertGreater(cosd, 0.8)
        self.assertEqual(o.select_comps(src), [3])

    def test_virac_grade(self):
        comps = []
        tgt = []
        for j in range(121):
            t = 55400.2 + 20 * j
            off = 0.1 if j % 2 else 0.0
            for sid in range(1, 7):
                comps.append({"sourceid": sid, "catid": 1000 + j, "t": t, "mag": 12.0 + 0.1 * sid + off})
            ecl = 2.0 if j in (13, 14, 15, 16) else 0.0
            tgt.append({"sourceid": o.VIRAC_ID, "catid": 1000 + j, "t": t, "mag": 12.5 + off + ecl})
        res, nights = o.virac_grade(tgt, comps)
        self.assertEqual(res["n_nights"], 121)
        self.assertEqual(res["n_baseline_nights"], 107)
        self.assertEqual(res["n_epochs_without_zp"], 0)
        self.assertAlmostEqual(res["baseline_ks"], 12.5, places=6)
        self.assertEqual(res["events"]["-1"]["n_nights_inside"], 4)
        self.assertEqual(res["events"]["-1"]["label"], "K_DEEP")
        self.assertAlmostEqual(res["events"]["-1"]["max_dm_inside"], 2.0, places=6)
        self.assertEqual(res["events"]["0"]["n_nights_inside"], 4)
        self.assertEqual(res["events"]["0"]["label"], "K_FLAT")
        self.assertEqual(res["label"], "K_DEEP")
        self.assertEqual(res["unexpected_dim_nights"], [])


def build(dim1=True, drop=(), extra_dim=()):
    target, tracer = [], []
    for i in range(66):
        t = 56875.3 + 30 * i
        for dt in (0.0, 0.07):
            tt = t + dt
            tracer.append((tt, 5000.0))
            if i in drop:
                continue
            f = 1000.0
            if i in (2, 3, 4) or i in extra_dim:
                f = 100.0
            if dim1 and i in (45, 46, 47):
                f = 100.0
            target.append((tt, f))
    return target, tracer


class TestGaia(unittest.TestCase):
    def test_parse_rows(self):
        rows = [
            {"source_id": "7", "band": "G", "time": "1000.5", "flux": "1500.0", "rejected_by_photometry": "False"},
            {"source_id": "7", "band": "BP", "time": "1000.5", "flux": "10"},
            {"source_id": "7", "band": "G", "time": "1001.5", "flux": "", "mag": "17.6874",
             "rejected_by_photometry": "false"},
            {"source_id": "7", "band": "G", "time": "1002.5", "flux": "1400", "rejected_by_photometry": "true"},
            {"source_id": "8", "band": "G", "time": "1003.5", "flux": "1400"},
        ]
        tr = o.parse_gaia_rows(rows, 7)
        self.assertEqual(len(tr), 2)
        self.assertAlmostEqual(tr[0][0], 56197.5)
        self.assertAlmostEqual(tr[0][1], 1500.0)
        self.assertAlmostEqual(tr[1][1], 10 ** 3.2, places=3)
        arr = [{"source_id": 7, "g_transit_time": "(1000.5, 1001.5, nan)",
                "g_transit_flux": "(1500.0, 1600.0, 1700.0)"}]
        tr2 = o.parse_gaia_rows(arr, 7)
        self.assertEqual(len(tr2), 2)
        self.assertAlmostEqual(tr2[1][0], 56198.5)
        self.assertAlmostEqual(tr2[1][1], 1600.0)

    def test_returned_on_schedule(self):
        target, tracer = build(dim1=True)
        res, pts = o.gaia_grade(target, tracer)
        self.assertEqual(res["n_baseline_transits"], 112)
        self.assertAlmostEqual(res["baseline_flux"], 1000.0)
        self.assertEqual(res["control"], "PASS")
        self.assertEqual(res["main"], "RETURNED_ON_SCHEDULE")
        self.assertEqual(res["event_2018"]["n_dim_transits"], 6)
        self.assertEqual(res["event_2014"]["n_dim_transits"], 6)
        self.assertEqual(res["off_window_dim_episodes"], [])
        p = res["period_interval_days"]
        self.assertAlmostEqual(p[0], 1260.07, places=2)
        self.assertAlmostEqual(p[1], 1319.93, places=2)
        self.assertEqual(len(res["forecast"]), 3)
        self.assertEqual(res["n_tracer_visits_missed"], 0)

    def test_not_returned(self):
        target, tracer = build(dim1=False)
        res, _ = o.gaia_grade(target, tracer)
        self.assertEqual(res["control"], "PASS")
        self.assertEqual(res["main"], "NOT_RETURNED")
        self.assertEqual(res["event_2018"]["n_bright_core_visits"], 2)

    def test_dropout_only(self):
        target, tracer = build(dim1=False, drop=(45, 46))
        res, _ = o.gaia_grade(target, tracer)
        self.assertEqual(res["event_2018"]["n_missed_core_visits"], 2)
        self.assertEqual(res["event_2018"]["n_bright_core_visits"], 0)
        self.assertEqual(res["main"], "DROPOUT_ONLY")

    def test_off_window_and_shifted(self):
        target, tracer = build(dim1=False, extra_dim=(10, 11))
        res, _ = o.gaia_grade(target, tracer)
        self.assertEqual(len(res["off_window_dim_episodes"]), 1)
        self.assertEqual(res["off_window_dim_episodes"][0]["n_dim"], 4)
        self.assertEqual(res["main"], "NOT_RETURNED")
        target, tracer = build(dim1=False, extra_dim=(50, 51))
        res, _ = o.gaia_grade(target, tracer)
        self.assertEqual(res["main"], "RETURNED_SHIFTED")
        self.assertAlmostEqual(res["shift_days"], 177.3, places=3)

    def test_control_fail(self):
        target, tracer = build(dim1=False)
        target = [(t, 1000.0) for t, f in target]
        res, _ = o.gaia_grade(target, tracer)
        self.assertEqual(res["control"], "FAIL")


if __name__ == "__main__":
    unittest.main()

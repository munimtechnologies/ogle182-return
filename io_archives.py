# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""Archive I/O for the OGLE-BLG182.1.162852 test (ESO TAP for VIRAC2; Gaia DataLink for DR4).

I/O only; the frozen rules live in ogle182.py. astropy is imported lazily (VOTable fallback only).
Made with AI agents working in MV2's name (MV2, Munim, Inc.).
"""
import csv
import glob
import io
import os
import time
import urllib.parse
import urllib.request
import zipfile

UA = {"User-Agent": "ogle182-return/0.1 (open research code; MV2, Munim, Inc.)"}
ESO_TAP = "https://archive.eso.org/tap_cat/sync"
GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
DATALINK = "https://gea.esac.esa.int/data-server/data?"


def fetch(url, timeout=600, tries=4, log=print):
    last = None
    for a in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as e:  # retry with back-off
            last = e
            log("fetch failed (%d/%d) %s: %s" % (a + 1, tries, url[:160], e))
            time.sleep(10 * (a + 1))
    raise RuntimeError("fetch failed: %s" % last)


def eso_tap_csv(query, maxrec=200000, log=print):
    q = urllib.parse.urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", MAXREC=str(maxrec),
                                    QUERY=query))
    text = fetch(ESO_TAP + "?" + q, log=log).decode("utf-8", "replace")
    if text.lstrip().startswith("<"):
        raise RuntimeError("ESO TAP returned XML (error?): %s" % text[:500])
    return list(csv.DictReader(io.StringIO(text)))


def virac_sources_box(ra, dec, half_deg_dec, cosd, log=print):
    dra = half_deg_dec / cosd
    q = ("SELECT sourceid, ra, de, phot_ks_mean_mag, ks_n_det, ks_mad, duplicate "
         "FROM safcat.VVVX_VIRAC_V2_SOURCES WHERE ra BETWEEN %.6f AND %.6f AND de BETWEEN %.6f AND %.6f"
         % (ra - dra, ra + dra, dec - half_deg_dec, dec + half_deg_dec))
    return eso_tap_csv(q, log=log)


def virac_lightcurves(source_ids, log=print):
    ids = ",".join(str(int(s)) for s in source_ids)
    q = ("SELECT sourceid, catid, mjdobs, filter, mag, emag, phot_flag, ambiguous_match, chi, seeing "
         "FROM safcat.VVVX_VIRAC_V2_LC WHERE sourceid IN (%s) AND filter = 'Ks'" % ids)
    return eso_tap_csv(q, log=log)


# --- Gaia (DR4 day) ---------------------------------------------------------------------------
def gaia_tap_csv(query, log=print):
    q = urllib.parse.urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=query))
    text = fetch(GAIA_TAP + "?" + q, log=log, tries=2).decode("utf-8", "replace")
    return list(csv.DictReader(io.StringIO(text)))


def datalink_urls(source_id, release="Gaia DR4"):
    out = []
    for rt in ("EPOCH_PHOTOMETRY", "EPOCH_PHOTOMETRY_G", "ALL"):
        for structure in ("DATAMODEL_STANDARD", "INDIVIDUAL", None):
            for ident in ("%s %d" % (release, source_id), str(source_id)):
                p = dict(RETRIEVAL_TYPE=rt, ID=ident, FORMAT="CSV", RELEASE=release)
                if structure:
                    p["DATA_STRUCTURE"] = structure
                out.append(DATALINK + urllib.parse.urlencode(p))
    return out


def read_table(path):
    ext = path.lower().rsplit(".", 1)[-1]
    if ext == "csv":
        with open(path, newline="", encoding="utf-8", errors="replace") as f:
            lines = [ln for ln in f if not ln.startswith("#")]
        return list(csv.DictReader(io.StringIO("".join(lines))))
    if ext in ("xml", "vot", "votable"):
        from astropy.io.votable import parse
        tab = parse(path).get_first_table().to_table(use_names_over_ids=True)
    else:
        from astropy.table import Table
        tab = Table.read(path)
    rows = []
    for r in tab:
        d = {}
        for n in tab.colnames:
            v = r[n]
            try:
                import numpy as np
                if np.ndim(v) > 0:
                    v = [None if m else x for x, m in zip(np.ma.asarray(v).data.tolist(),
                                                          np.ma.getmaskarray(np.ma.asarray(v)).tolist())]
                elif isinstance(v, (bytes, np.bytes_)):
                    v = v.decode()
                elif isinstance(v, np.generic):
                    v = v.item()
            except ImportError:
                pass
            d[n] = v
        rows.append(d)
    return rows


def save_response(blob, dest_dir, stem, log=print):
    """Store a DataLink response (zip, csv or votable); return readable table paths."""
    os.makedirs(dest_dir, exist_ok=True)
    head = blob[:200].lstrip().lower()
    if blob[:2] == b"PK":
        zp = os.path.join(dest_dir, stem + ".zip")
        with open(zp, "wb") as f:
            f.write(blob)
        out = []
        try:
            with zipfile.ZipFile(zp) as z:
                for name in z.namelist():
                    if name.endswith("/"):
                        continue
                    tgt = os.path.join(dest_dir, stem + "__" + os.path.basename(name))
                    with open(tgt, "wb") as f:
                        f.write(z.read(name))
                    out.append(tgt)
        except zipfile.BadZipFile as e:
            log("bad zip %s: %s" % (zp, e))
        return out
    if head.startswith(b"<!doctype html") or head.startswith(b"<html") or not blob.strip():
        log("rejected %s: HTML or empty" % stem)
        return []
    ext = "xml" if head.startswith(b"<?xml") or head.startswith(b"<votable") else "csv"
    p = os.path.join(dest_dir, stem + "." + ext)
    with open(p, "wb") as f:
        f.write(blob)
    return [p]


def gaia_epoch_rows(source_id, dest_dir, release="Gaia DR4", log=print):
    """Try the DataLink variants until one gives rows; return (rows, url) or ([], None)."""
    for i, url in enumerate(datalink_urls(source_id, release)):
        try:
            blob = fetch(url, tries=2, log=log)
        except Exception as e:
            log("datalink failed: %s" % e)
            continue
        for p in save_response(blob, dest_dir, "%d_try%02d" % (source_id, i), log=log):
            try:
                rows = read_table(p)
            except Exception as e:
                log("unreadable %s: %s" % (p, str(e)[:200]))
                continue
            if rows and ("g_transit_time" in rows[0] or "band" in rows[0]):
                log("epoch photometry: %d rows from %s" % (len(rows), url))
                return rows, url
            log("not epoch photometry: %s columns %s" % (p, sorted(rows[0].keys())[:12] if rows else []))
    return [], None


def cached_files(d):
    return sorted(glob.glob(os.path.join(d, "*")))

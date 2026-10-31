#!/usr/bin/env python3
"""Build server/cuzk_dmp1g_index.json.gz: the map sheets of ČÚZK's DMP 1G (the Czech Republic's digital surface model,
1st generation: airborne lidar 2009-2013, LAZ in S-JTSK / EPSG:5514, CC BY 4.0) for the Czech bridge decks
(server/roadz.py "cz", BNE_CZ_DECKS).

The source is ČÚZK's ATOM download feed (~24 MB, 16,301 sheets in 2026-10). Its sheets are the SM5 layout: 2,500 x
2,000 m cells of S-JTSK; each entry's WGS84 box gives its cell (the centre, through sjtsk.forward), kept as
"col,row" -> sheet name (col = floor(x / 2500), row = floor(y / 2000)). The sheet's file is
openzu.cuzk.cz/opendata/DMP1G/epsg-5514/<name>.zip (one LAZ, ~5-10 MB).

    python3 tools/build_cz_dmp_index.py DMP1G-SJTSK.xml     # a downloaded copy
    python3 tools/build_cz_dmp_index.py                     # download it first
"""
import gzip
import json
import os
import re
import sys
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "server"))
import sjtsk  # noqa: E402

URL = "https://atom.cuzk.cz/DMP1G-SJTSK/DMP1G-SJTSK.xml"
OUT = os.path.join(ROOT, "server", "cuzk_dmp1g_index.json.gz")
SHEET_W, SHEET_H = 2500.0, 2000.0
UA = "BeamNG.Earth-build (+https://github.com/TecnaGamer/BeamNG.Earth)"


def main():
    if len(sys.argv) > 1:
        text = open(sys.argv[1], encoding="utf-8").read()
    else:
        req = urllib.request.Request(URL, headers={"User-Agent": UA})
        text = urllib.request.urlopen(req, timeout=300).read().decode("utf-8")
    tiles = {}
    for e in text.split("<entry")[1:]:
        n = re.search(r"DMP1G-SJTSK_([A-Z0-9]+)\.xml", e)
        p = re.search(r"<georss:polygon>([^<]*)", e)
        if not n or not p:
            continue
        v = [float(t) for t in p.group(1).split()]
        la, lo = v[0::2], v[1::2]
        x, y = sjtsk.forward((min(la) + max(la)) / 2, (min(lo) + max(lo)) / 2)
        k = f"{int(x // SHEET_W)},{int(y // SHEET_H)}"
        if k in tiles:
            raise SystemExit(f"{n.group(1)} and {tiles[k]} in one cell {k}")
        tiles[k] = n.group(1)
    with gzip.open(OUT, "wt", encoding="utf-8") as fh:
        json.dump({"source": URL, "sheet": [SHEET_W, SHEET_H],
                   "url": "https://openzu.cuzk.cz/opendata/DMP1G/epsg-5514/{name}.zip",
                   "tiles": dict(sorted(tiles.items()))}, fh, separators=(",", ":"))
    print(f"{len(tiles)} sheets -> {OUT} ({os.path.getsize(OUT) // 1024} KB)")


if __name__ == "__main__":
    main()

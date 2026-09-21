#!/usr/bin/env python3
"""
Step 11a (atlas): stream ImmuNexUT nominal eQTL (E-GEAD-420; 43 GB zip holding a tar.gz)
and keep only rows whose variant lies in the atlas windows (hg38) from Step 10.
One TSV per cell type in OUT_DIR; resumable (cells listed in _done.txt are skipped).

usage:  python 11a_immunexut_stream_atlas.py  URL  PROJECT_DIR
        URL = the same E-GEAD-420 zip URL used for the first extraction
"""
import sys, os, time, tarfile, warnings
warnings.filterwarnings("ignore")
url, proj = sys.argv[1], os.path.expanduser(sys.argv[2])
win_f = os.path.join(proj, "results", "atlas", "windows_hg38.tsv")
out_dir = os.path.join(proj, "data", "ImmuNexUT", "atlas_hg38"); os.makedirs(out_dir, exist_ok=True)

BY_CHR = {}
with open(win_f) as fh:
    next(fh)                                           # header: name chr start end
    for line in fh:
        name, c, s, e = line.rstrip("\n").split("\t")[:4]
        c = c if c.startswith("chr") else "chr" + c
        BY_CHR.setdefault(c.encode(), []).append((int(s), int(e), name.encode()))
print(f"{sum(len(v) for v in BY_CHR.values())} windows on {len(BY_CHR)} chromosomes", flush=True)

def open_member(u):
    if os.path.exists(u):
        import zipfile; z = zipfile.ZipFile(u)
    else:
        from remotezip import RemoteZip
        z = RemoteZip(u, initial_buffer_size=20 * 1024 * 1024, timeout=600)
    m = z.infolist()[0]
    return z.open(m.filename)

done_path = os.path.join(out_dir, "_done.txt")
done = set(open(done_path).read().split()) if os.path.exists(done_path) else set()
t0 = time.time(); nmem = 0
with tarfile.open(fileobj=open_member(url), mode="r|gz") as tar:
    for ti in tar:
        if not ti.isfile():
            continue
        nmem += 1
        cell = os.path.basename(ti.name).replace("_nominal.txt", "")
        if cell in done:
            print(f"[{nmem}] {cell}: already done, skipping", flush=True); continue
        fh = tar.extractfile(ti); header = fh.readline()
        tmp = os.path.join(out_dir, cell + ".tsv.part"); kept = nread = 0; t1 = time.time()
        with open(tmp, "wb") as out:
            out.write(b"window\t" + header)
            for line in fh:
                nread += 1
                p = line.split(b"\t", 11)
                w = BY_CHR.get(p[8])
                if w is not None:
                    pos = int(p[9])
                    for s, e, name in w:
                        if s <= pos <= e:
                            out.write(name + b"\t" + line); kept += 1; break
                if nread % 20_000_000 == 0:
                    print(f"   {cell}: {nread/1e6:.0f}M rows, kept {kept}, {(time.time()-t0)/60:.0f} min", flush=True)
        os.replace(tmp, os.path.join(out_dir, cell + ".tsv"))
        with open(done_path, "a") as d: d.write(cell + "\n")
        print(f"[{nmem}] {cell}: {nread/1e6:.1f}M rows -> kept {kept} ({(time.time()-t1)/60:.1f} min; total {(time.time()-t0)/60:.0f} min)", flush=True)
print(f"ALL DONE: {nmem} cell types in {(time.time()-t0)/60:.0f} min", flush=True)

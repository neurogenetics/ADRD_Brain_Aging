import os
import csv
import concurrent.futures
import synapseclient
import synapseutils

syn = synapseclient.login()  # or syn.login(authToken="...")

PARENT = "syn64425506"  # folder to crawl
DEST = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/seaad/src_data"  # where to save chosen files
OUT_MAP = f"{DEST}/file_synid_map.tsv"  # filename -> synID table


# Extensions you actually want. Tune as needed.
WANTED_EXT = (".tsv.gz", ".tsv.gz.tbi", ".h5")

# 1. Crawl the folder (recursively) and build filename -> synID map.
#    walk() yields (dirpath_tuple, subdirs, files) where files is a
#    list of (filename, synID) tuples.
rows = []
for _dirpath, _dirnames, files in synapseutils.walk(syn, PARENT, includeTypes=["file"]):
    for fname, synid in files:
        rows.append((fname, synid))

# 2. Save the full mapping so you can inspect / subset it however you like.
with open(OUT_MAP, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["filename", "synid"])
    w.writerows(rows)
print(f"Wrote {len(rows)} entries to {OUT_MAP}")

# 3. Subset to just the file types you want.
wanted = [(name, synid) for name, synid in rows if name.lower().endswith(WANTED_EXT)]
print(f"{len(wanted)} files match {WANTED_EXT}")


# 4. Download only the subset.
def download_file(file_info):
    name, synid = file_info
    print(f"Starting download for {name} ({synid})...")
    try:
        syn.get(synid, downloadLocation=DEST)
        print(f"Finished downloading {name} ({synid})")
    except Exception as e:
        print(f"Error downloading {name} ({synid}): {e}")


print(f"Beginning parallel download of {len(wanted)} files...")
# You can increase max_workers based on network capabilities
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
    executor.map(download_file, wanted)
print("All downloads completed.")

# 5. Verify and retry missing files (1 attempt).
missing = []
for name, synid in wanted:
    # Synapse usually retains the original filename
    file_path = os.path.join(DEST, name)
    if not os.path.exists(file_path):
        missing.append((name, synid))

if missing:
    print(f"Found {len(missing)} missing files. Retrying downloads...")
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
        executor.map(download_file, missing)

    # Final check
    still_missing = []
    for name, synid in missing:
        if not os.path.exists(os.path.join(DEST, name)):
            still_missing.append(name)

    if still_missing:
        print(f"Failed to download {len(still_missing)} files after retry:")
        for m in still_missing:
            print(f"  - {m}")
    else:
        print("All files successfully verified after retry!")
else:
    print("All files successfully verified on the first pass!")

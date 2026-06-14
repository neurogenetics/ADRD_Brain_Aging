# https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyDetails?Study=syn26223298
# https://www.synapse.org/Synapse:syn64425506
#instructions for download in README; access token in laptop documents (created interactively on synapse website)
cd /data/ADRD/brain_aging/exploration/data/Phase2/SEAAD_EC/
module load synapseclient

###to only get out the things I want
synapse login
#will prompt to put in access token
# synapse config # only need to do this the first time I'm downloading. Do if `synapse login` doesn't work

# synapse get -r syn64425506 #this would also download all their processed data. which I don't really want.

## the code below downloads every fastq file from this project
# python <<'EOF' > fastq_manifest.tsv
# import synapseclient
# import synapseutils
# 
# syn = synapseclient.login(silent=True)
# 
# for _, _, files in synapseutils.walk(
#         syn,
#         "syn64425506",
#         includeTypes=["file"]):
#     for fname, synid in files:
#         if fname.lower().endswith(
#             (".fastq.gz", ".fastq", ".fq.gz", ".fq")
#         ):
#             print(f"{synid}\t{fname}")
# EOF
# 
# wc -l fastq_manifest.tsv
# 
# mkdir -p SEAAD_FASTQ
# 
# cut -f1 fastq_manifest.tsv | while read synid
# do
#     synapse get "$synid" --downloadLocation SEAAD_FASTQ
# done

###this allegedly has to be done interactively and looks like it will take ~5 hr to download all (even though biowulf is going strong at ~140MB/s, maybe I'm missing somethign that says I should use helix?). If I cannot do in one sitting I need to pause the process and check which are downloaded, subset fastq_manifest down and continue download.

##### need to model what is below this line to make fastq_manifest_remaining1.tsv - based on what we already downloaded and what to do next

# fastq_manifest_remaining1.tsv #<- how to check files we have and see what is missing from fastq_manifest.tsv
# #NY-AT0648_7 was the last one that was starting, it does I1, R1, R2, R3 successively. the atac for 17 samples downloaded in 5 hours. next need to do the other 17-23. followed by the RNA. 2.187 Tb into what I think approximate close to 5-8?

echo "Files: $(find . -type f | wc -l)"
du -BG -d 1 .

# I'm going to make this fastq_manifest_remaining1 file in an interactive job and try to proceed in helix.

# ## the code below downloads every fastq file from this project
# python <<'EOF'
# from pathlib import Path
# 
# manifest_path = Path("fastq_manifest.tsv")
# download_dir = Path("SEAAD_FASTQ")
# out_path = Path("fastq_manifest_remaining1.tsv")
# 
# # Collect downloaded FASTQ basenames from the local download directory
# downloaded = set()
# if download_dir.exists():
#     for p in download_dir.rglob("*"):
#         if p.is_file() and p.name.lower().endswith((".fastq.gz", ".fastq", ".fq.gz", ".fq")):
#             downloaded.add(p.name)
# 
# remaining = []
# with manifest_path.open() as fh:
#     for line in fh:
#         line = line.rstrip("\n")
#         if not line:
#             continue
#         synid, fname = line.split("\t", 1)
#         if fname not in downloaded:
#             remaining.append((synid, fname))
# 
# with out_path.open("w") as oh:
#     for synid, fname in remaining:
#         oh.write(f"{synid}\t{fname}\n")
# 
# print(f"Downloaded FASTQ files found locally: {len(downloaded)}")
# print(f"Manifest entries total: {sum(1 for _ in manifest_path.open())}")
# print(f"Remaining entries written: {len(remaining)}")
# print(f"Wrote: {out_path}")
# EOF
# 
# wc -l fastq_manifest_remaining1.tsv
# #run this on helix?
# cut -f1 fastq_manifest_remaining1.tsv | while read -r synid
# do
#     synapse get "$synid" --downloadLocation SEAAD_FASTQ
# done

#ugghh still takes forever, also helix is not much faster...
#try to swarm by prefix 
WORKDIR="/data/ADRD/brain_aging/exploration/data/Phase2/SEAAD_EC"
SCRIPT_DIR="/data/ADRD/brain_aging/exploration/scripts/Phase2/AD"
MANIFEST="${WORKDIR}/fastq_manifest_remaining1.tsv"
DOWNLOAD_DIR="${WORKDIR}/SEAAD_FASTQ"
GROUP_DIR="${WORKDIR}/group_ids"
SWARM_FILE="${SCRIPT_DIR}/2_download_groups.swarm"

mkdir -p "${DOWNLOAD_DIR}" "${GROUP_DIR}"

python <<'EOF'
from pathlib import Path
from collections import defaultdict
import re

workdir = Path("/data/ADRD/brain_aging/exploration/data/Phase2/SEAAD_EC")
group_dir = workdir / "group_ids"
manifest = workdir / "fastq_manifest_remaining1.tsv"
swarm_file = Path("/data/ADRD/brain_aging/exploration/scripts/Phase2/AD/2_download_groups.swarm")

groups = defaultdict(list)

with manifest.open() as fh:
    for line in fh:
        line = line.rstrip("\n")
        if not line:
            continue
        synid, fname = line.split("\t", 1)

        base = Path(fname).name
        prefix = re.sub(
            r'_(?:I\d+|R\d+)_[0-9]{3}\.fastq(?:\.gz)?$',
            '',
            base,
            flags=re.IGNORECASE
        )

        groups[prefix].append((synid, fname))

for prefix in sorted(groups):
    out = group_dir / f"{prefix}.tsv"
    with out.open("w") as oh:
        for synid, fname in groups[prefix]:
            oh.write(f"{synid}\t{fname}\n")

with swarm_file.open("w") as sh:
    for prefix in sorted(groups):
        group_tsv = group_dir / f"{prefix}.tsv"
        sh.write(
            "bash -lc 'set -euo pipefail; "
            "cd /data/ADRD/brain_aging/exploration/data/Phase2/SEAAD_EC; "
            f"cut -f1 {group_tsv} | while read -r synid; do "
            "synapse get \"$synid\" --downloadLocation /data/ADRD/brain_aging/exploration/data/Phase2/SEAAD_EC/SEAAD_FASTQ; "
            "done'\n"
        )

print(f"Wrote {len(groups)} groups")
print(f"Wrote swarm file: {swarm_file}")
EOF

wc -l "${SWARM_FILE}"
head -n 3 "${SWARM_FILE}"

swarm -f /data/ADRD/brain_aging/exploration/scripts/Phase2/AD/2_download_groups.swarm -g 4 -t 1 --module synapseclient

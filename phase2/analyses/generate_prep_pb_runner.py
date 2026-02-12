import stat
from pathlib import Path

# Constants
QUANTS_DIR = Path("/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/quants")
PROJECT_PREFIX = "aging_phase2"


def generate_job_script():
    if not QUANTS_DIR.exists():
        print(f"Error: Quants directory not found at {QUANTS_DIR}")
        return

    # Find all parquet files matching the pattern
    files = list(QUANTS_DIR.glob(f"{PROJECT_PREFIX}.*.*.parquet"))

    prep_commands_by_modality = {"rna": [], "atac": []}
    vp_commands_by_modality = {"rna": [], "atac": []}

    print(f"Scanning {QUANTS_DIR} for input files...")

    for file_path in files:
        # Expected format: aging_phase2.{cell_type}.{modality}.parquet
        parts = file_path.name.split(".")

        if len(parts) < 4:
            print(f"Skipping malformed filename: {file_path.name}")
            continue

        modality = parts[-2]

        # Validate modality
        if modality not in ["rna", "atac"]:
            continue

        cell_type_parts = parts[1:-2]
        cell_type = ".".join(cell_type_parts)

        # Prep PB Data Command
        prep_cmd = (
            f"uv run phase2/analyses/prep_pb_data.py "
            f"--modality {modality} "
            f"--cell-type {cell_type}"
        )
        prep_commands_by_modality[modality].append(prep_cmd)

        # Variance Partition Command
        vp_cmd = (
            f"uv run phase2/analyses/run_variance_partition.py "
            f"--modality {modality} "
            f"--cell-type {cell_type}"
        )
        vp_commands_by_modality[modality].append(vp_cmd)

    # Function to write and chmod script
    def write_script(output_path, commands, description):
        if not commands:
            return
        
        commands.sort()
        with open(output_path, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(f"# Generated script to {description}\n\n")

            for cmd in commands:
                f.write(f'echo "Running: {cmd}"\n')
                f.write(f"{cmd}\n")
                f.write("\n")
        
        st = output_path.stat()
        output_path.chmod(st.st_mode | stat.S_IEXEC)
        print(f"Successfully generated {output_path} with {len(commands)} commands.")

    # Write prep scripts
    for modality, commands in prep_commands_by_modality.items():
        write_script(
            Path(f"run_prep_pb_jobs_{modality}.sh"),
            commands,
            f"run prep_pb_data.py for {modality} modality"
        )

    # Write variance partition scripts
    for modality, commands in vp_commands_by_modality.items():
        write_script(
            Path(f"run_variance_partition_jobs_{modality}.sh"),
            commands,
            f"run run_variance_partition.py for {modality} modality"
        )


if __name__ == "__main__":
    generate_job_script()

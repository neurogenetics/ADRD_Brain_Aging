import stat
from pathlib import Path

# Constants
QUANTS_DIR = Path("/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/quants")
OUTPUT_SCRIPT = Path("run_prep_pb_jobs.sh")
PROJECT_PREFIX = "aging_phase2"

def generate_job_script():
    if not QUANTS_DIR.exists():
        print(f"Error: Quants directory not found at {QUANTS_DIR}")
        return

    # Find all parquet files matching the pattern
    files = list(QUANTS_DIR.glob(f"{PROJECT_PREFIX}.*.*.parquet"))
    
    commands_by_modality = {"rna": [], "atac": []}
    
    print(f"Scanning {QUANTS_DIR} for input files...")
    
    for file_path in files:
        # Expected format: aging_phase2.{cell_type}.{modality}.parquet
        parts = file_path.name.split('.')
        
        if len(parts) < 4:
            print(f"Skipping malformed filename: {file_path.name}")
            continue
            
        modality = parts[-2]
        
        # Validate modality
        if modality not in ["rna", "atac"]:
            continue
            
        cell_type_parts = parts[1:-2]
        cell_type = ".".join(cell_type_parts)
        
        cmd = (
            f"uv run phase2/analyses/prep_pb_data.py "
            f"--modality {modality} "
            f"--cell-type {cell_type}"
        )
        commands_by_modality[modality].append(cmd)

    # Write separate scripts for each modality
    for modality, commands in commands_by_modality.items():
        if not commands:
            continue
            
        commands.sort()
        output_script = Path(f"run_prep_pb_jobs_{modality}.sh")
        
        with open(output_script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(f"# Generated script to run prep_pb_data.py for {modality} modality\n\n")
            
            # Add basic error handling or logging to the bash script
            f.write("set -e\n\n")
            
            for cmd in commands:
                f.write(f'echo "Running: {cmd}"\n')
                f.write(f"{cmd}\n")
                f.write("\n")

        # Make the script executable
        st = output_script.stat()
        output_script.chmod(st.st_mode | stat.S_IEXEC)

        print(f"Successfully generated {output_script} with {len(commands)} commands.")

if __name__ == "__main__":
    generate_job_script()

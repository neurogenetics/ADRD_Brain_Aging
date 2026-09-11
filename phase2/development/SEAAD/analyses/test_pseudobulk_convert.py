import os
import sys
import tempfile
import unittest
import numpy as np
import pandas as pd
import anndata as ad
import mudata as md
from subprocess import run


class TestDiseasePseudobulkConvert(unittest.TestCase):
    def setUp(self):
        # Create temp directory
        self.test_dir = tempfile.TemporaryDirectory()

        # Create 120 cells: 40 cells per donor (donor_1, donor_2, donor_3)
        # For each donor: 20 Neurons, 20 Astrocytes
        donors = ["donor_1"] * 40 + ["donor_2"] * 40 + ["donor_3"] * 40
        cell_types = (["Neuron"] * 20 + ["Astrocyte"] * 20) * 3
        obs_index = [f"cell_{i}" for i in range(120)]

        # Create dummy AnnData objects with some count data
        rna_data = ad.AnnData(
            X=np.random.randint(10, 100, size=(120, 10)).astype(np.float32),
            obs=pd.DataFrame({
                "sample_id": donors,
                "broad_cell_type": cell_types
            }, index=obs_index),
            var=pd.DataFrame(index=[f"gene_{i}" for i in range(10)])
        )

        atac_data = ad.AnnData(
            X=np.random.randint(0, 10, size=(120, 15)).astype(np.float32),
            obs=pd.DataFrame({
                "sample_id": donors,
                "broad_cell_type": cell_types
            }, index=obs_index),
            var=pd.DataFrame(index=[f"peak_{i}" for i in range(15)])
        )

        # Create MuData object
        self.mdata = md.MuData({"rna": rna_data, "atac": atac_data})
        self.mdata_path = os.path.join(self.test_dir.name, "test_input.h5mu")
        self.mdata.write(self.mdata_path)

        self.quants_dir = os.path.join(self.test_dir.name, "quants")
    def tearDown(self):
        self.test_dir.cleanup()
        
    def test_cli_pseudobulk_convert(self):
        # Define command line arguments
        cmd = [
            sys.executable,
            "phase2/development/SEAAD/analyses/pseudobulk_convert.py",
            "-i", self.mdata_path,
            "-o", self.quants_dir,
            "--project", "test_project",
            "--work-dir", self.test_dir.name,
            "--cell-type-col", "broad_cell_type",
            "--sample-col", "sample_id",
            "--aggregate-type", "sum",
        ]
        
        # Run CLI script
        result = run(cmd, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, f"Script failed with output:\n{result.stderr}\n{result.stdout}")
        
        # Verify output parquet files exist
        # Output format: {project_name}.{cell_type}.{modality}.parquet
        # Cell types: Neuron, Astrocyte
        expected_files = [
            "test_project.Neuron.rna.parquet",
            "test_project.Neuron.atac.parquet",
            "test_project.Astrocyte.rna.parquet",
            "test_project.Astrocyte.atac.parquet",
        ]
        
        for fname in expected_files:
            file_path = os.path.join(self.quants_dir, fname)
            self.assertTrue(os.path.exists(file_path), f"Expected file not found: {file_path}")
            
            # Read and verify content
            df = pd.read_parquet(file_path)
            self.assertGreater(df.shape[0], 0, f"Empty dataframe in {fname}")
            self.assertGreater(df.shape[1], 0, f"No columns in {fname}")
            
            # Rows should be sample IDs (index: donor_1, donor_2, etc.)
            self.assertIn("donor_1", df.index)
            
        # Verify log file exists
        log_path = os.path.join(self.test_dir.name, "logs", "test_project_disease_pseudobulk_convert.log")
        self.assertTrue(os.path.exists(log_path))
        with open(log_path, "r") as f:
            log_content = f.read()
        self.assertIn("Starting Disease Pseudobulk Conversion", log_content)
        self.assertIn("Saved", log_content)
        self.assertIn("Disease Pseudobulk Conversion Finished Successfully", log_content)


if __name__ == "__main__":
    unittest.main()

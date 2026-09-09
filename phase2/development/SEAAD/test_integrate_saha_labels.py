import os
import sys
import tempfile
import unittest
import numpy as np
import pandas as pd
import anndata as ad
import mudata as md
from subprocess import run


class TestIntegrateSahaLabels(unittest.TestCase):
    def setUp(self):
        # Create temp directory
        self.test_dir = tempfile.TemporaryDirectory()
        
        # Create dummy AnnData objects
        obs_index = ["cell_A", "cell_B", "cell_C"]
        rna_data = ad.AnnData(
            X=np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]),
            obs=pd.DataFrame(index=obs_index),
            var=pd.DataFrame(index=["gene_1", "gene_2"])
        )
        atac_data = ad.AnnData(
            X=np.array([[0.0, 1.0], [1.0, 0.0], [0.0, 0.0]]),
            obs=pd.DataFrame(index=obs_index),
            var=pd.DataFrame(index=["peak_1", "peak_2"])
        )
        
        # Create MuData object
        self.mdata = md.MuData({"rna": rna_data, "atac": atac_data})
        self.mdata_path = os.path.join(self.test_dir.name, "test_input.h5mu")
        self.mdata.write(self.mdata_path)
        
        # Create dummy cell labels CSV
        # One cell has different labels, one has duplicate entries (to test deduplication)
        self.labels_path = os.path.join(self.test_dir.name, "test_labels.csv")
        labels_data = {
            "barcode": ["cell_A", "cell_B", "cell_C", "cell_B"],  # cell_B is duplicated
            "cell_type": ["Neuron", "Astrocyte", "Microglia", "Astrocyte"],
            "subtype": ["Excitatory", "SST", "PVALB", "SST"],
        }
        pd.DataFrame(labels_data).to_csv(self.labels_path, index=False)
        
        self.output_path = os.path.join(self.test_dir.name, "test_output.h5mu")
        
    def tearDown(self):
        self.test_dir.cleanup()
        
    def test_cli_integration(self):
        # Define command line arguments
        cmd = [
            sys.executable,
            "phase2/development/SEAAD/integrate_saha_labels.py",
            "-i", self.mdata_path,
            "-o", self.output_path,
            "-l", self.labels_path,
            "-b", "barcode",
            "-c", "cell_type", "subtype",
        ]
        
        # Run CLI script
        result = run(cmd, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, f"Script failed with output:\n{result.stderr}\n{result.stdout}")
        
        # Verify output exists
        self.assertTrue(os.path.exists(self.output_path))
        
        # Read and check integrated data
        mdata_out = md.read(self.output_path)
        
        # Check global .obs
        self.assertIn("cell_type", mdata_out.obs.columns)
        self.assertIn("subtype", mdata_out.obs.columns)
        
        self.assertEqual(mdata_out.obs.loc["cell_A", "cell_type"], "Neuron")
        self.assertEqual(mdata_out.obs.loc["cell_B", "cell_type"], "Astrocyte")
        self.assertEqual(mdata_out.obs.loc["cell_C", "cell_type"], "Microglia")
        
        self.assertEqual(mdata_out.obs.loc["cell_A", "subtype"], "Excitatory")
        self.assertEqual(mdata_out.obs.loc["cell_B", "subtype"], "SST")
        self.assertEqual(mdata_out.obs.loc["cell_C", "subtype"], "PVALB")
        
        # Verify log file exists and is populated with correct information
        log_path = os.path.join(self.test_dir.name, "test_output.log")
        self.assertTrue(os.path.exists(log_path))
        with open(log_path, "r") as f:
            log_content = f.read()
        self.assertIn("Logging to terminal and file", log_content)
        self.assertIn("Alignment Statistics:", log_content)
        self.assertIn("Cells in MuData updated with a label: 3 (100.00%)", log_content)
        self.assertIn("Labels from CSV matched and added to MuData: 3 (100.00%)", log_content)
        self.assertIn("Value counts for categorical column 'cell_type':", log_content)
        self.assertIn("Neuron", log_content)
        self.assertIn("Astrocyte", log_content)
        self.assertIn("Microglia", log_content)

        # Check that individual modalities also got the columns mapped correctly
        for mod_name in ["rna", "atac"]:
            mod_obs = mdata_out.mod[mod_name].obs
            self.assertIn("cell_type", mod_obs.columns)
            self.assertIn("subtype", mod_obs.columns)
            self.assertEqual(mod_obs.loc["cell_A", "cell_type"], "Neuron")
            self.assertEqual(mod_obs.loc["cell_B", "subtype"], "SST")


if __name__ == "__main__":
    unittest.main()

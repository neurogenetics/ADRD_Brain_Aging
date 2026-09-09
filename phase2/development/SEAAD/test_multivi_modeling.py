import os
import sys
import tempfile
import unittest
import numpy as np
import pandas as pd
import anndata as ad
import mudata as md
from subprocess import run


class TestRunMultiVIModeling(unittest.TestCase):
    def setUp(self):
        # Create temp directory
        self.test_dir = tempfile.TemporaryDirectory()
        
        # Create dummy AnnData objects
        obs_index = [f"cell_{i}" for i in range(100)]
        
        # Ensure var names are unique and appropriate
        rna_var_names = [f"gene_{i}" for i in range(20)]
        atac_var_names = [f"peak_{i}" for i in range(30)]
        
        rna_data = ad.AnnData(
            X=np.random.randint(0, 10, size=(100, 20)).astype(np.float32),
            obs=pd.DataFrame(index=obs_index),
            var=pd.DataFrame(index=rna_var_names)
        )
        rna_data.layers["counts"] = rna_data.X.copy()
        
        atac_data = ad.AnnData(
            X=np.random.randint(0, 2, size=(100, 30)).astype(np.float32),
            obs=pd.DataFrame(index=obs_index),
            var=pd.DataFrame(index=atac_var_names)
        )
        atac_data.layers["counts"] = atac_data.X.copy()
        
        # Create MuData object
        self.mdata = md.MuData({"rna": rna_data, "atac": atac_data})
        
        # Add metadata fields to global obs
        self.mdata.obs["cell_type"] = np.random.choice(["Neuron", "Astrocyte"], size=100)
        self.mdata.obs["batch_id"] = np.random.choice(["batch_1", "batch_2"], size=100)
        
        self.mdata_path = os.path.join(self.test_dir.name, "test_input.h5mu")
        self.mdata.write(self.mdata_path)
        
        self.output_path = os.path.join(self.test_dir.name, "test_output.h5mu")
        
    def tearDown(self):
        self.test_dir.cleanup()
        
    def test_cli_modeling(self):
        # Define command line arguments
        cmd = [
            sys.executable,
            "phase2/development/SEAAD/multivi_modeling.py",
            "-i", self.mdata_path,
            "-o", self.output_path,
            "-c", "cell_type", "batch_id",
            "-b", "batch_id",
            "-e", "2",      # Only train 2 epochs for fast test execution
            "-l", "3",      # Small latent dimension
            "--detect-hv-features",
            "--top-genes", "10",
            "--top-peaks", "15",
        ]
        
        # Run CLI script
        result = run(cmd, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, f"Script failed with output:\n{result.stderr}\n{result.stdout}")
        
        # Verify output mudata exists
        self.assertTrue(os.path.exists(self.output_path))
        
        # Read and check modeled data
        mdata_out = md.read(self.output_path)
        
        # Verify that features were subsetted correctly (smaller than original shapes)
        self.assertLess(mdata_out.n_vars, 50)
        self.assertLess(mdata_out.mod["rna"].n_vars, 20)
        self.assertLess(mdata_out.mod["atac"].n_vars, 30)
        
        # Check that MultiVI and UMAP embeddings were created
        self.assertIn("X_multivi", mdata_out.obsm.keys())
        self.assertIn("X_umap", mdata_out.obsm.keys())
        
        # Verify shape of coordinates matches expectations (100 cells, 3 dimensions for multivi latent)
        self.assertEqual(mdata_out.obsm["X_multivi"].shape, (100, 3))
        self.assertEqual(mdata_out.obsm["X_umap"].shape, (100, 2))
        
        # Verify log file exists and is populated
        log_path = os.path.join(self.test_dir.name, "test_output.log")
        self.assertTrue(os.path.exists(log_path))
        with open(log_path, "r") as f:
            log_content = f.read()
        
        # If test is going to fail, print log content first for debugging
        plot_celltype_path = os.path.join(self.test_dir.name, "test_output_umap_cell_type.png")
        if not os.path.exists(plot_celltype_path):
            print("--- LOG FILE CONTENT ---")
            print(log_content)
            print("------------------------")

        self.assertIn("Starting MultiVI Modeling Pipeline", log_content)
        self.assertIn("Found 'counts' layer in RNA modality, will use for MultiVI model setup.", log_content)
        self.assertIn("Found 'counts' layer in ATAC modality, will use for MultiVI model setup.", log_content)
        self.assertIn("Training MultiVI model for up to 2 epochs", log_content)
        self.assertIn("MultiVI Modeling Pipeline Finished Successfully", log_content)
        
        # Verify that UMAP figure plots were saved
        plot_celltype_path = os.path.join(self.test_dir.name, "test_output_umap_cell_type.png")
        plot_batch_path = os.path.join(self.test_dir.name, "test_output_umap_batch_id.png")
        
        self.assertTrue(os.path.exists(plot_celltype_path))
        self.assertTrue(os.path.exists(plot_batch_path))
        
        # Ensure figure files are populated (not empty)
        self.assertGreater(os.path.getsize(plot_celltype_path), 0)
        self.assertGreater(os.path.getsize(plot_batch_path), 0)


if __name__ == "__main__":
    unittest.main()

import os
import sys
import tempfile
import unittest
import numpy as np
import pandas as pd
from subprocess import run


class TestDiseaseRegressionPipeline(unittest.TestCase):
    def setUp(self):
        # Create temp directory
        self.test_dir = tempfile.TemporaryDirectory()
        
        self.info_dir = os.path.join(self.test_dir.name, "sample_info")
        self.quants_dir = os.path.join(self.test_dir.name, "quants")
        self.results_dir = os.path.join(self.test_dir.name, "results")
        self.figs_dir = os.path.join(self.test_dir.name, "figures")
        self.logs_dir = os.path.join(self.test_dir.name, "logs")
        
        os.makedirs(self.info_dir, exist_ok=True)
        os.makedirs(self.quants_dir, exist_ok=True)
        os.makedirs(self.results_dir, exist_ok=True)
        os.makedirs(self.figs_dir, exist_ok=True)
        os.makedirs(self.logs_dir, exist_ok=True)
        
        # We need mock final_covariates for two cell types: Neuron and Astrocyte
        self.cell_types = ["Neuron", "Astrocyte"]
        self.modalities = ["rna", "atac"]
        self.project = "seaad_ec_multiome"
        
        # Create 10 mock donors (donor_1 through donor_10)
        donors = [f"donor_{i}" for i in range(1, 11)]
        
        for ct in self.cell_types:
            for modality in self.modalities:
                # 1. Mock final covariates (contains 'dx', 'cell_counts', sex, race, etc. + PCA)
                covars_data = {
                    "dx": [0, 0, 0, 0, 0, 1, 1, 1, 1, 1], # 5 Control, 5 AD
                    "sex": ["M", "F", "M", "F", "M", "F", "M", "F", "M", "F"],
                    "race": ["Caucasian"] * 10,
                    "ageDeath": [75.0, 82.5, 91.0, 78.0, 88.0, 84.0, 92.5, 79.0, 85.0, 90.0],
                    "PMI": [5.5, 6.2, 4.8, 8.0, 7.2, 5.0, 9.5, 6.0, 7.5, 8.2],
                    "pH": [6.5, 6.8, 6.2, 6.6, 6.4, 6.3, 6.7, 6.5, 6.6, 6.4],
                    "brainWeight": [1250.0, 1310.0, 1180.0, 1420.0, 1290.0, 1200.0, 1350.0, 1260.0, 1310.0, 1280.0],
                    "cell_counts": [100, 150, 80, 200, 120, 90, 140, 110, 160, 130],
                    f"{ct}_counts": [100, 150, 80, 200, 120, 90, 140, 110, 160, 130],
                    "PCA_0": np.random.randn(10).tolist(),
                    "PCA_1": np.random.randn(10).tolist(),
                }
                
                cov_df = pd.DataFrame(covars_data, index=donors)
                cov_file = os.path.join(self.info_dir, f"{self.project}.{ct}.{modality}.final_covariates.csv")
                cov_df.to_csv(cov_file)
                
                # 2. Mock quants parquets
                n_features = 10 if modality == "rna" else 15
                cols = [f"gene_{i}" if modality == "rna" else f"peak_{i}" for i in range(n_features)]
                
                # Make some features slightly correlated with dx to ensure we get nominal significance
                quants_matrix = np.random.rand(10, n_features).astype(np.float32)
                # Introduce a direct disease association in the first feature
                quants_matrix[5:, 0] += 2.0
                
                quant_df = pd.DataFrame(quants_matrix, index=donors, columns=cols)
                quant_file = os.path.join(self.quants_dir, f"{self.project}.{ct}.{modality}.parquet")
                quant_df.to_parquet(quant_file)

    def tearDown(self):
        self.test_dir.cleanup()
        
    def test_disease_regression_and_downstream_aux(self):
        # 1. Run pseudobulk_regression.py for both cell types (Neuron & Astrocyte)
        for ct in self.cell_types:
            cmd_reg = [
                sys.executable,
                "phase2/development/SEAAD/analyses/pseudobulk_regression.py",
                "--project", self.project,
                "--work-dir", self.test_dir.name,
                "--modality", "rna",
                "--cell-type", ct,
                "--target-variable", "dx",
                "--regression-type", "ols",
            ]
            res_reg = run(cmd_reg, capture_output=True, text=True)
            self.assertEqual(res_reg.returncode, 0, f"pseudobulk_regression failed for {ct}:\n{res_reg.stderr}\n{res_reg.stdout}")
            
            # Also run a robust RLM (vwrlm) regression to support the filter_regression_type_differences script
            cmd_robust = [
                sys.executable,
                "phase2/development/SEAAD/analyses/pseudobulk_regression.py",
                "--project", self.project,
                "--work-dir", self.test_dir.name,
                "--modality", "rna",
                "--cell-type", ct,
                "--target-variable", "dx",
                "--regression-type", "vwrlm",
            ]
            res_robust = run(cmd_robust, capture_output=True, text=True)
            self.assertEqual(res_robust.returncode, 0, f"pseudobulk_regression robust failed for {ct}:\n{res_robust.stderr}\n{res_robust.stdout}")

        # Assert output OLS and VWRLM files exist for both Neuron and Astrocyte
        for ct in self.cell_types:
            self.assertTrue(os.path.exists(os.path.join(self.results_dir, f"{self.project}.rna.{ct}.ols.dx.csv")))
            self.assertTrue(os.path.exists(os.path.join(self.results_dir, f"{self.project}.rna.{ct}.vwrlm.dx.csv")))

        # 2. Run post_pseudobulk_regression.py for both types (ols and vwrlm)
        for rtype in ["ols", "vwrlm"]:
            cmd_post = [
                sys.executable,
                "phase2/development/SEAAD/analyses/post_pseudobulk_regression.py",
                "--project", self.project,
                "--work-dir", self.test_dir.name,
                "--modality", "rna",
                "--target-variable", "dx",
                "--regression-type", rtype,
            ]
            res_post = run(cmd_post, capture_output=True, text=True)
            self.assertEqual(res_post.returncode, 0, f"post_pseudobulk_regression failed for {rtype}:\n{res_post.stderr}\n{res_post.stdout}")

        # Assert post-processed outputs exist
        ols_full_file = os.path.join(self.results_dir, f"{self.project}.all_celltypes.rna.ols.dx.csv")
        ols_fdr_file = os.path.join(self.results_dir, f"{self.project}.all_celltypes.rna.ols_fdr.dx.csv")
        self.assertTrue(os.path.exists(ols_full_file))
        self.assertTrue(os.path.exists(ols_fdr_file))
        self.assertTrue(os.path.exists(os.path.join(self.results_dir, f"{self.project}.all_celltypes.rna.vwrlm.dx.csv")))

        # Overwrite the FDR results file with mock significant rows to guarantee downstream power analysis works
        mock_sig_fdr_df = pd.DataFrame({
            "feature": ["gene_0", "gene_1", "gene_2", "gene_0", "gene_1", "gene_2"],
            "tissue": ["Neuron", "Neuron", "Neuron", "Astrocyte", "Astrocyte", "Astrocyte"],
            "coef": [1.2, -0.8, 0.5, 1.1, -0.9, 0.6],
            "p-value": [0.001, 0.002, 0.005, 0.0015, 0.003, 0.004],
            "fdr_bh": [0.01, 0.02, 0.03, 0.015, 0.025, 0.035],
            "percentchange": [120.0, -80.0, 50.0, 110.0, -90.0, 60.0],
            "fc": [2.2, 0.2, 1.5, 2.1, 0.1, 1.6],
            "log2fc": [1.1, -1.2, 0.5, 1.0, -1.3, 0.6],
            "stderr": [0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
            "z": [12.0, -8.0, 5.0, 11.0, -9.0, 6.0],
        })
        mock_sig_fdr_df.to_csv(ols_fdr_file, index=False)

        # 3. Run filter_regression_type_differences.py
        cmd_filter = [
            sys.executable,
            "phase2/development/SEAAD/analyses/filter_regression_type_differences.py",
            "--project", self.project,
            "--work-dir", self.test_dir.name,
            "--modality", "rna",
            "--target-variable", "dx",
            "--general-type", "ols",
            "--robust-type", "vwrlm",
        ]
        res_filter = run(cmd_filter, capture_output=True, text=True)
        self.assertEqual(res_filter.returncode, 0, f"filter_regression_type_differences failed:\n{res_filter.stderr}\n{res_filter.stdout}")
        
        # Assert clean filtered output was written
        filtered_results_file = os.path.join(self.results_dir, f"{self.project}.rna.all_celltypes.ols_fdr_filtered.dx.csv")
        self.assertTrue(os.path.exists(filtered_results_file))
        
        # Overwrite the filtered results with mock significant rows to guarantee downstream plot generation works
        mock_sig_df = pd.DataFrame({
            "feature": ["gene_0", "gene_1", "gene_2", "gene_0", "gene_1", "gene_2"],
            "tissue": ["Neuron", "Neuron", "Neuron", "Astrocyte", "Astrocyte", "Astrocyte"],
            "coef": [1.2, -0.8, 0.5, 1.1, -0.9, 0.6],
            "p-value": [0.001, 0.002, 0.005, 0.0015, 0.003, 0.004],
            "fdr_bh": [0.01, 0.02, 0.03, 0.015, 0.025, 0.035],
            "percentchange": [120.0, -80.0, 50.0, 110.0, -90.0, 60.0],
            "fc": [2.2, 0.2, 1.5, 2.1, 0.1, 1.6],
            "log2fc": [1.1, -1.2, 0.5, 1.0, -1.3, 0.6],
            "stderr": [0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
            "z": [12.0, -8.0, 5.0, 11.0, -9.0, 6.0],
        })
        mock_sig_df.to_csv(filtered_results_file, index=False)

        # 4. Run run_variance_partition.py for 'Neuron' cell type
        cmd_var = [
            sys.executable,
            "phase2/development/SEAAD/analyses/run_variance_partition.py",
            "--project", self.project,
            "--work-dir", self.test_dir.name,
            "--modality", "rna",
            "--cell-type", "Neuron",
            "--target-variable", "dx",
        ]
        res_var = run(cmd_var, capture_output=True, text=True)
        self.assertEqual(res_var.returncode, 0, f"run_variance_partition failed:\n{res_var.stderr}\n{res_var.stdout}")

        # Assert variance boxplots & partitions exist
        self.assertTrue(os.path.exists(os.path.join(self.figs_dir, f"{self.project}_Neuron_rna_variance_partition_known.csv")))
        self.assertTrue(os.path.exists(os.path.join(self.figs_dir, f"{self.project}_Neuron_rna_variance_boxen_known.png")))

        # 5. Run cell_counts_regression.py
        cmd_counts = [
            sys.executable,
            "phase2/development/SEAAD/analyses/cell_counts_regression.py",
            "--project", self.project,
            "--work-dir", self.test_dir.name,
            "--target-variable", "dx",
            "--weight-term", "cell_counts",
        ]
        res_counts = run(cmd_counts, capture_output=True, text=True)
        self.assertEqual(res_counts.returncode, 0, f"cell_counts_regression failed:\n{res_counts.stderr}\n{res_counts.stdout}")
        self.assertTrue(os.path.exists(os.path.join(self.results_dir, f"{self.project}.all_cell_counts_bias.dx.csv")))

        # 6. Run celltype_target_effect_similarity.py (figures)
        cmd_similarity = [
            sys.executable,
            "phase2/development/SEAAD/figures/celltype_target_effect_similarity.py",
            "--project", self.project,
            "--work-dir", self.test_dir.name,
            "--modality", "rna",
            "--target-variable", "dx",
            "--regression-type", "ols",
        ]
        res_similarity = run(cmd_similarity, capture_output=True, text=True)
        self.assertEqual(res_similarity.returncode, 0, f"celltype_target_effect_similarity failed:\n{res_similarity.stderr}\n{res_similarity.stdout}")
        
        # Assert clustermap image exists
        plot_path = os.path.join(self.figs_dir, f"{self.project}_rna_dx_effect_similarity_clustermap.png")
        if not os.path.exists(plot_path):
            print("--- SIMILARITY SCRIPT STDOUT ---")
            print(res_similarity.stdout)
            print("--- SIMILARITY SCRIPT STDERR ---")
            print(res_similarity.stderr)
            print("--------------------------------")
        self.assertTrue(os.path.exists(plot_path))

        # 7. Run regression_power_analysis.py
        # Mock some results to look like they have significant values for power curves
        # Let's specify the paths explicitly using --results
        cmd_power = [
            sys.executable,
            "phase2/development/SEAAD/analyses/regression_power_analysis.py",
            "--project", self.project,
            "--work-dir", self.test_dir.name,
            "--target-variable", "dx",
            "--regression-type", "ols",
            "--labels", "RNA",
            "--sizes", "10",
            "--results", ols_fdr_file,
            "--output", os.path.join(self.figs_dir, "dx_WLS_Power_Curve.png"),
        ]
        res_power = run(cmd_power, capture_output=True, text=True)
        self.assertEqual(res_power.returncode, 0, f"regression_power_analysis failed:\n{res_power.stderr}\n{res_power.stdout}")
        
        # Assert power curves saved successfully
        self.assertTrue(os.path.exists(os.path.join(self.figs_dir, "dx_WLS_Power_Curve.png")))


if __name__ == "__main__":
    unittest.main()

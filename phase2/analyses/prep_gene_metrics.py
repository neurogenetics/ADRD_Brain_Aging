import argparse
import logging
from pathlib import Path
import sys
import csv

import pandas as pd
import numpy as np

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(description="Build gene metrics from UCSC knownGene CSV.")
    parser.add_argument(
        "--input", 
        type=str, 
        default="/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/ucsc/hg38_knownGene.csv",
        help="Path to the input UCSC knownGene CSV file."
    )
    parser.add_argument(
        "--output", 
        type=str, 
        default="/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/ucsc/hg38_gene_metrics.csv",
        help="Path to save the output metrics CSV."
    )
    return parser.parse_args()

def main():
    args = parse_args()
    input_file = Path(args.input)
    output_file = Path(args.output)
    
    if not input_file.exists():
        logger.error(f"Input file not found: {input_file}")
        sys.exit(1)
        
    logger.info(f"Loading data from {input_file}...")
    
    # Custom parser because exonStarts/exonEnds contain comma-separated arrays which break standard CSV parsers
    parsed_data = []
    with open(input_file, 'r', encoding='utf-8') as f:
        # We don't use csv.reader directly because the internal commas aren't quoted. 
        # But looking at UCSC format, the columns before exonStarts are fixed (8 cols), 
        # and the columns after exonEnds are fixed (5 cols).
        # We can extract based on exonCount.
        
        header = f.readline().strip().split(',')
        # 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'proteinID', 'alignID', 'kgID', 'geneSymbol', 'refseq'
        
        for i, line in enumerate(f):
            parts = line.strip().split(',')
            
            try:
                name = parts[0]
                txStart = int(parts[3])
                txEnd = int(parts[4])
                exonCount = int(parts[7])
                
                # Exon arrays have `exonCount` number of elements
                starts_end_idx = 8 + exonCount
                exonStarts = [int(x) for x in parts[8:starts_end_idx] if x]
                
                ends_end_idx = starts_end_idx + exonCount
                # +1 handles the trailing comma in the exonEnds list
                exonEnds = [int(x) for x in parts[starts_end_idx:ends_end_idx+1] if x]
                
                # The trailing columns might have empty strings at the end if the last few columns are blank
                # We count backwards from the end of the line. 
                # According to header: ..., 'kgID', 'geneSymbol', 'refseq'
                geneSymbol = parts[-2]
                
                # Skip if geneSymbol is empty
                if not geneSymbol:
                    continue
                
                if len(exonStarts) != exonCount or len(exonEnds) != exonCount:
                    continue
                    
                # Calculate tx_gene_size
                tx_gene_size = sum(end - start for start, end in zip(exonStarts, exonEnds))
                tx_genomic_region_size = txEnd - txStart
                
                parsed_data.append({
                    "name": name,
                    "geneSymbol": geneSymbol,
                    "exonCount": exonCount,
                    "tx_gene_size": tx_gene_size,
                    "tx_genomic_region_size": tx_genomic_region_size
                })
            except Exception as e:
                # Silently skip malformed lines
                continue

    df = pd.DataFrame(parsed_data)
    logger.info(f"Successfully parsed {len(df)} transcripts.")

    logger.info("Aggregating metrics per geneSymbol...")
    
    # 1. Number of transcripts: unique values in 'name'
    transcript_counts = df.groupby("geneSymbol")["name"].nunique().rename("num_transcripts")
    
    # 2. Number of exons: max 'exonCount'
    max_exons = df.groupby("geneSymbol")["exonCount"].max().rename("max_exons")
    
    # 3. Gene size: max 'tx_gene_size'
    max_gene_size = df.groupby("geneSymbol")["tx_gene_size"].max().rename("max_gene_size")
    
    # 4. Genomic region size: max 'tx_genomic_region_size'
    max_genomic_size = df.groupby("geneSymbol")["tx_genomic_region_size"].max().rename("max_genomic_region_size")
    
    # Combine into a single DataFrame
    metrics_df = pd.concat([
        transcript_counts,
        max_exons,
        max_gene_size,
        max_genomic_size
    ], axis=1)
    
    # Ensure index is named for feature_enrichment.py
    metrics_df.index.name = "feature"
    
    logger.info(f"Generated metrics for {len(metrics_df)} unique genes.")
    
    output_file.parent.mkdir(parents=True, exist_ok=True)
    metrics_df.to_csv(output_file)
    logger.info(f"Saved gene metrics to {output_file}")

if __name__ == "__main__":
    main()
extract_features_from_gb.py

This script generates a phylogeny-ready locus set from chloroplast GenBank files. It extracts coding and non-coding features (CDS, tRNA, rRNA, exons, introns, and intergenic spacers). Outputs are organised in staged directories: 1_with_copies (raw extractions), 2_without_copies (filtered), 3_combined (per-locus multi-FASTAs), and a final 1_for_phylo_analysis set with standardised numeric prefixes (1_ CDS/exon, 2_ rRNA, 3_ tRNA, 4_ intron, 5_ spacer). A summary.csv file will be generated, recording feature type, name, coordinates, strand, length, and source file.

Usage: python extract_features_from_gb.py <input_gb_dir> <output_dir>

Recommended: Python 3.9+ and Biopython ≥1.84.

pip install biopython

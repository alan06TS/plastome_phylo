**extract_features_from_gb.py**

This script generates a phylogeny-ready locus set from chloroplast GenBank files. It extracts coding and non-coding features (CDS, tRNA, rRNA, exons, introns, and intergenic spacers). Outputs are organised in staged directories: 1_with_copies (raw extractions), 2_without_copies (filtered), 3_combined (per-locus multi-FASTAs), and a final 1_for_phylo_analysis set with standardised numeric prefixes (1_ CDS/exon, 2_ rRNA, 3_ tRNA, 4_ intron, 5_ spacer). A summary.csv file will be generated, recording feature type, name, coordinates, strand, length, and source file.

Usage: python /path/to/extract_features_from_gb.py <input_gb_dir_path> <output_dir_path>

Recommended: Python 3.9+ and Biopython ≥1.84.

pip install biopython

Citation: Alan T.S., Xiang C. L., Krishnaraj T.P. and Sunojkumar P. (2025) Characterisation and phylogenetic analysis of the complete chloroplast genome of the rare endemic Isodon rivularis (Lamiaceae) from the Western Ghats, Journal of Asia-Pacific Biodiversity, (in press).

**IMPORTANT NOTE**: Check the generated FASTA files for errors. Check transplicing genes, for example, the rps12 gene. Check start and stop codons of protein-coding genes. Remember to avoid the use of synonyms of gene names, for example, rrn4.5, sometimes rrn4.5S; hence, use rrn4.5 in all GenBank files.

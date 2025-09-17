#!/usr/bin/env python3
import os
import csv
import sys
import shutil
from pathlib import Path
from collections import defaultdict

from Bio import SeqIO, SeqFeature
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def create_output_dirs(base_dir):
    dirs = {}
    for key in [
        "cds", "exons", "introns", "trna", "rrna",
        "intergenic_spacers", "pseudogenes"
    ]:
        path = base_dir / key / "1_with_copies"
        path.mkdir(parents=True, exist_ok=True)
        dirs[key] = path
    return dirs

def clean_filename(name: str) -> str:
    """Replace '-' and '.' with '_' in filenames (excluding extensions)."""
    return name.replace("-", "_").replace(".", "_")

def write_fasta(seq, label, feature_type, parent_dir, gb_stem):
    cleaned_label = clean_filename(label)
    subfolder = parent_dir / gb_stem
    subfolder.mkdir(parents=True, exist_ok=True)
    filepath = subfolder / f"{cleaned_label}.fasta"
    record = SeqRecord(seq, id=gb_stem, description="")
    with open(filepath, "w") as f:
        SeqIO.write(record, f, "fasta")

def find_nearest_gene(position, gene_coords, direction):
    """Find nearest upstream/downstream gene from a list of (start,end,name,strand)."""
    if direction == 'upstream':
        candidates = [g for g in gene_coords if g[1] <= position]
        return max(candidates, key=lambda x: x[1], default=None)
    else:  # downstream
        candidates = [g for g in gene_coords if g[0] >= position]
        return min(candidates, key=lambda x: x[0], default=None)

def extract_intergenic(seq, start, end, label, gb_stem,
                       output_dirs, summary_rows,
                       source_file, extracted_starts):
    """Extract and write an intergenic spacer only if its start coordinate hasn’t already been used."""
    if end <= start or start in extracted_starts:
        return
    spacer_seq = seq[start:end]
    write_fasta(spacer_seq, label, "intergenic_spacers",
                output_dirs["intergenic_spacers"], gb_stem)
    summary_rows.append([
        "intergenic_spacer",
        label,
        start,
        end,
        "+",
        len(spacer_seq),
        source_file
    ])
    extracted_starts.add(start)

def concatenate_exons(exon_dir, gb_stem):
    exon_folder = exon_dir / gb_stem
    if not exon_folder.exists():
        return
    exon_files = list(exon_folder.glob("*_exon*.fasta"))
    exons_by_gene = defaultdict(list)
    for f in exon_files:
        parts = f.stem.split("_exon")
        if len(parts) == 2:
            gene = parts[0]
            try:
                exon_num = int(parts[1])
            except ValueError:
                continue
            exons_by_gene[gene].append((exon_num, f))

    for gene, exon_list in exons_by_gene.items():
        full_seq = Seq("")
        for _, f in sorted(exon_list, key=lambda x: x[0]):
            record = next(SeqIO.parse(f, "fasta"))
            full_seq += record.seq
        final_record = SeqRecord(full_seq, id=gb_stem, description="")
        output_path = exon_folder / f"{clean_filename(gene)}.fasta"
        with open(output_path, "w") as out:
            SeqIO.write(final_record, out, "fasta")

def copy_filtered_files(output_dir):
    """
    Copy from 1_with_copies to 2_without_copies (or 2_pseudo for pseudogenes).
    """
    for category in [
        "cds", "exons", "introns", "trna", "rrna",
        "intergenic_spacers", "pseudogenes"
    ]:
        src = output_dir / category / "1_with_copies"
        # Special case: pseudogenes → 2_pseudo
        if category == "pseudogenes":
            dst = output_dir / category / "2_pseudo"
        else:
            dst = output_dir / category / "2_without_copies"

        dst.mkdir(parents=True, exist_ok=True)
        log_path = dst / "log.txt"
        counter = 1
        with open(log_path, "w") as log:
            for gb_folder in src.glob("*"):
                if not gb_folder.is_dir():
                    continue

                start_spacer_path = None
                end_spacer_path = None
                if category == "intergenic_spacers":
                    for fasta in gb_folder.glob("*.fasta"):
                        stem_lc = fasta.stem.lower()
                        if stem_lc == "start_to_first_gene_spacer":
                            start_spacer_path = fasta
                        elif stem_lc == "last_gene_to_end_spacer":
                            end_spacer_path = fasta

                for fasta in gb_folder.glob("*.fasta"):
                    name = fasta.stem
                    name_lc = name.lower()

                    if category == "pseudogenes":
                        # just copy everything, skip filtering
                        pass
                    elif category != "intergenic_spacers":
                        if "exon" in name or name_lc.count("copy") >= 1:
                            log.write(f"Skipped (non-intergenic filter): {fasta.name}\n")
                            continue
                    else:
                        if name_lc in {"start_to_first_gene_spacer", "last_gene_to_end_spacer"}:
                            log.write(f"Skipped (edge spacer, to concatenate): {fasta.name}\n")
                            continue
                        if name_lc.count("copy") >= 2:
                            log.write(f"Skipped (>=2 'copy' in name): {fasta.name}\n")
                            continue

                    subfolder = dst / name
                    subfolder.mkdir(parents=True, exist_ok=True)
                    new_name = f"{clean_filename(name)}_{counter}.fasta"
                    shutil.copy(fasta, subfolder / new_name)
                    log.write(f"Copied: {fasta.name} -> {subfolder}/{new_name}\n")
                    counter += 1

                if category == "intergenic_spacers":
                    if start_spacer_path is not None and end_spacer_path is not None:
                        rec_end = next(SeqIO.parse(end_spacer_path, "fasta"))
                        rec_start = next(SeqIO.parse(start_spacer_path, "fasta"))
                        combined_seq = rec_end.seq + rec_start.seq
                        combined_record = SeqRecord(combined_seq, id=gb_folder.name, description="")

                        subfolder = dst / "end_gene_to_start_gene_spacer"
                        subfolder.mkdir(parents=True, exist_ok=True)
                        new_name = f"{clean_filename('end_gene_to_start_gene_spacer')}_{counter}.fasta"
                        combined_path = subfolder / new_name
                        with open(combined_path, "w") as out:
                            SeqIO.write(combined_record, out, "fasta")
                        log.write(
                            f"Created concatenated spacer: {combined_path.name} "
                            f"from {end_spacer_path.name} + {start_spacer_path.name}\n"
                        )
                        counter += 1
                    else:
                        log.write(
                            f"Concatenation skipped for {gb_folder.name}: "
                            f"missing start_to_first_gene_spacer or last_gene_to_end_spacer.\n"
                        )

            if category == "intergenic_spacers":
                for subfolder in dst.glob("*"):
                    if not subfolder.is_dir():
                        continue
                    for fasta in subfolder.glob("*.fasta"):
                        if fasta.stem.lower().count("copy") >= 2:
                            fasta.unlink()
                            with open(log_path, "a") as log2:
                                log2.write(f"Deleted duplicate (>=2 'copy'): {fasta.name}\n")

def combine_to_3_combined(output_dir):
    for category in [
        "cds", "exons", "introns", "trna", "rrna",
        "intergenic_spacers", "pseudogenes"
    ]:
        # Adjust for pseudogenes → 2_pseudo
        if category == "pseudogenes":
            src = output_dir / category / "2_pseudo"
        else:
            src = output_dir / category / "2_without_copies"

        combined_path = output_dir / category / "3_combined"
        combined_path.mkdir(parents=True, exist_ok=True)
        if not src.exists():
            continue
        for subfolder in src.glob("*"):
            if not subfolder.is_dir():
                continue
            records = []
            for f in subfolder.glob("*.fasta"):
                records.extend(list(SeqIO.parse(f, "fasta")))
            if records:
                combined_file = combined_path / f"{clean_filename(subfolder.name)}.fasta"
                SeqIO.write(records, combined_file, "fasta")

def create_phylo_folder(output_dir):
    phylo_path = output_dir / "1_for_phylo_analysis"
    phylo_path.mkdir(parents=True, exist_ok=True)
    prefix_map = {
        "cds": "1_",
        "rrna": "2_",
        "trna": "3_",
        "introns": "4_",
        "intergenic_spacers": "5_"
        # pseudogenes intentionally excluded
    }

    # Copy for normal categories except exons
    for cat, prefix in prefix_map.items():
        src = output_dir / cat / "3_combined"
        if not src.exists():
            continue
        for f in src.glob("*.fasta"):
            new_name = f"{prefix}{clean_filename(f.stem)}.fasta"
            shutil.copy(f, phylo_path / new_name)

    # Special handling for exons: ensure ALL exons transferred
    exon_src = output_dir / "exons" / "3_combined"
    if exon_src.exists():
        for f in exon_src.glob("*.fasta"):
            base = f.stem
            base_lc = base.lower()
            if base_lc.startswith("trn"):
                prefix = "3_"
            elif base_lc.startswith("rrn"):
                prefix = "2_"
            else:
                prefix = "1_"
            new_name = f"{prefix}{clean_filename(base)}.fasta"
            shutil.copy(f, phylo_path / new_name)

def get_outer_gene_coords(gene_coords):
    """Keep only intervals not fully nested within another."""
    sorted_coords = sorted(gene_coords, key=lambda x: (x[0], -x[1]))
    outer = []
    for start, end, name, strand in sorted_coords:
        if not any(o_start <= start and o_end >= end for o_start, o_end, _, _ in outer):
            outer.append((start, end, name, strand))
    return outer

def process_gb_file(gb_file, output_dirs, summary_rows):
    record = SeqIO.read(gb_file, "genbank")
    sequence = record.seq
    gb_stem = gb_file.stem

    extracted_spacer_starts = set()

    # ---------------------------
    # 1) Build initial coords for nearest-gene lookups (include pseudogene-only 'gene' features)
    # Use display names WITHOUT the word 'pseudogene' for spacer labeling.
    # ---------------------------
    initial_gene_counts = defaultdict(int)          # for CDS/tRNA/rRNA copies
    initial_pseudo_label_counts = defaultdict(int)  # for pseudogene boundary label copies
    initial_gene_coords = []

    for feat in record.features:
        has_pseudo = "pseudo" in feat.qualifiers

        # Standard coding/RNA features
        if feat.type in ["CDS", "tRNA", "rRNA"]:
            gene = feat.qualifiers.get("gene", ["unknown"])[0]
            initial_gene_counts[gene] += 1
            name = (f"{gene}_copy{initial_gene_counts[gene]}" if initial_gene_counts[gene] > 1 else gene)
            s, e = int(feat.location.start), int(feat.location.end)
            initial_gene_coords.append((s, e, name, feat.location.strand))
            continue

        # If a /pseudo 'gene' feature exists (no CDS/tRNA/rRNA), include as boundary
        if has_pseudo and feat.type == "gene":
            gene = feat.qualifiers.get("gene", ["unknown"])[0]
            initial_pseudo_label_counts[gene] += 1
            # Label like other genes for spacers (NO 'pseudogene' in name)
            label_name = (f"{gene}_copy{initial_pseudo_label_counts[gene]}"
                          if initial_pseudo_label_counts[gene] > 1 else gene)
            s, e = int(feat.location.start), int(feat.location.end)
            initial_gene_coords.append((s, e, label_name, feat.location.strand))

    # ---------------------------
    # 2) Iterate features to extract sequences
    # ---------------------------
    gene_counts = defaultdict(int)                # CDS/tRNA/rRNA copy counts
    pseudogene_counts = defaultdict(int)          # for pseudogene FASTA naming
    pseudogene_label_counts = defaultdict(int)    # for spacer label naming (no 'pseudogene')
    gene_coords = []  # used to compute intergenic spacers across the whole molecule

    for feat in record.features:
        # handle pseudogenes (any feature with /pseudo)
        if "pseudo" in feat.qualifiers:
            gene = feat.qualifiers.get("gene", ["unknown"])[0]
            strand = feat.location.strand
            s, e = int(feat.location.start), int(feat.location.end)
            subseq = sequence[s:e]
            if strand == -1:
                subseq = subseq.reverse_complement()

            # FASTA name for pseudogene sequences: keep 'pseudogene'
            pseudogene_counts[gene] += 1
            pseudo_seq_name = f"{gene}_pseudogene_copy{pseudogene_counts[gene]}"
            write_fasta(subseq, pseudo_seq_name, "pseudogenes", output_dirs["pseudogenes"], gb_stem)
            summary_rows.append(["pseudogene", pseudo_seq_name, s, e, strand, len(subseq), gb_file.name])

            # Spacer boundary label for this pseudogene: like other genes (no 'pseudogene')
            pseudogene_label_counts[gene] += 1
            label_name = (f"{gene}_copy{pseudogene_label_counts[gene]}"
                          if pseudogene_label_counts[gene] > 1 else gene)
            gene_coords.append((s, e, label_name, strand))
            continue

        # skip anything not CDS/tRNA/rRNA
        if feat.type not in ["CDS", "tRNA", "rRNA"]:
            continue

        gene = feat.qualifiers.get("gene", ["unknown"])[0]
        gene_counts[gene] += 1
        name = (f"{gene}_copy{gene_counts[gene]}" if gene_counts[gene] > 1 else gene)
        strand = feat.location.strand
        ftype = feat.type.lower()

        # Special rps12 handling (unchanged)
        if gene == "rps12" and isinstance(feat.location, SeqFeature.CompoundLocation):
            parts = feat.location.parts
            exons = []
            for part in parts:
                ps, pe = int(part.start), int(part.end)
                seq = sequence[ps:pe]
                if part.strand == -1:
                    seq = seq.reverse_complement()
                exons.append({"start": ps, "end": pe,
                              "length": pe-ps,
                              "strand": part.strand,
                              "seq": seq})
            sorted_by_start = sorted(exons, key=lambda x: x["start"])
            first = sorted_by_start[0]
            largest = max(exons, key=lambda x: x["length"])
            remaining = [e for e in exons if e not in (first, largest)][0]
            ordered = [first, largest, remaining]

            for i, ex in enumerate(ordered, 1):
                lbl = f"{name}_exon{i}"
                write_fasta(ex["seq"], lbl, "exons", output_dirs["exons"], gb_stem)
                summary_rows.append([
                    "rps12_exon", lbl,
                    ex["start"], ex["end"],
                    ex["strand"], ex["length"], gb_file.name
                ])

            ex2, ex3 = ordered[1], ordered[2]
            s, e = sorted([ex2["end"], ex3["start"]])
            intr_seq = sequence[s:e]
            if strand == -1:
                intr_seq = intr_seq.reverse_complement()
            intr_lbl = f"{name}_intron1"
            write_fasta(intr_seq, intr_lbl, "introns", output_dirs["introns"], gb_stem)
            summary_rows.append([
                "intron", intr_lbl,
                s, e, strand,
                len(intr_seq), gb_file.name
            ])

            for ex in (ordered[0], ordered[2]):
                up = find_nearest_gene(ex["start"], initial_gene_coords, "upstream")
                if up and up[1] < ex["start"]:
                    extract_intergenic(
                        sequence, up[1], ex["start"],
                        f"{up[2]}_to_{name}_exon_spacer",
                        gb_stem, output_dirs, summary_rows,
                        gb_file.name, extracted_spacer_starts
                    )
                down = find_nearest_gene(ex["end"], initial_gene_coords, "downstream")
                if down and ex["end"] < down[0]:
                    extract_intergenic(
                        sequence, ex["end"], down[0],
                        f"{name}_exon_to_{down[2]}_spacer",
                        gb_stem, output_dirs, summary_rows,
                        gb_file.name, extracted_spacer_starts
                    )
            continue

        # Default feature writing for CDS/tRNA/rRNA
        s, e = int(feat.location.start), int(feat.location.end)
        subseq = sequence[s:e]
        if strand == -1:
            subseq = subseq.reverse_complement()
        write_fasta(subseq, name, ftype, output_dirs[ftype], gb_stem)
        summary_rows.append([ftype, name, s, e, strand, len(subseq), gb_file.name])
        gene_coords.append((s, e, name, strand))

        # If multi-exon, emit exons and introns
        if ftype in ["cds", "trna"] and isinstance(feat.location, SeqFeature.CompoundLocation):
            parts = sorted(feat.location.parts,
                           key=lambda x: int(x.start),
                           reverse=(strand == -1))
            for i, part in enumerate(parts, 1):
                ps, pe = int(part.start), int(part.end)
                exseq = sequence[ps:pe]
                if strand == -1:
                    exseq = exseq.reverse_complement()
                lbl = f"{name}_exon{i}"
                write_fasta(exseq, lbl, "exons", output_dirs["exons"], gb_stem)
                summary_rows.append(["exon", lbl, ps, pe, strand, len(exseq), gb_file.name])
            for i in range(len(parts)-1):
                a_end = int(parts[i].end)
                b_start = int(parts[i+1].start)
                if strand == -1:
                    a_end, b_start = int(parts[i+1].end), int(parts[i].start)
                if b_start > a_end:
                    intrseq = sequence[a_end:b_start]
                    if strand == -1:
                        intrseq = intrseq.reverse_complement()
                    lbl = f"{name}_intron{i+1}"
                    write_fasta(intrseq, lbl, "introns", output_dirs["introns"], gb_stem)
                    summary_rows.append([
                        "intron", lbl,
                        a_end, b_start, strand,
                        len(intrseq), gb_file.name
                    ])

    # ---------------------------
    # 3) Compute intergenic spacers using ALL boundaries (incl. pseudogenes),
    #    but spacer labels never include the word 'pseudogene'
    # ---------------------------
    outer_coords = get_outer_gene_coords(gene_coords)
    for i in range(len(outer_coords)-1):
        s1, e1, n1, _ = outer_coords[i]
        s2, _, n2, _ = outer_coords[i+1]
        extract_intergenic(
            sequence, e1, s2,
            f"{n1}_to_{n2}",
            gb_stem, output_dirs, summary_rows,
            gb_file.name, extracted_spacer_starts
        )
    if outer_coords:
        first_s = outer_coords[0][0]
        last_e  = outer_coords[-1][1]
        extract_intergenic(
            sequence, 0, first_s,
            "start_to_first_gene_spacer",
            gb_stem, output_dirs, summary_rows,
            gb_file.name, extracted_spacer_starts
        )
        extract_intergenic(
            sequence, last_e, len(sequence),
            "last_gene_to_end_spacer",
            gb_stem, output_dirs, summary_rows,
            gb_file.name, extracted_spacer_starts
        )

    # concatenate any exon pieces written earlier
    concatenate_exons(output_dirs["exons"], gb_stem)

def main():
    if len(sys.argv) != 3:
        print("Usage: python extract_features_from_gb.py <input_gb_dir> <output_dir>")
        sys.exit(1)

    input_dir = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])
    if not input_dir.exists() or not input_dir.is_dir():
        print(f"Error: Input directory '{input_dir}' is not valid.")
        sys.exit(1)

    output_dir.mkdir(exist_ok=True, parents=True)
    output_dirs = create_output_dirs(output_dir)
    summary_rows = []

    for gb_file in input_dir.glob("*.gb"):
        print(f"🔍 Processing: {gb_file.name}")
        process_gb_file(gb_file, output_dirs, summary_rows)

    # global summary
    with open(output_dir / "summary.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Feature_Type", "Gene_Name",
            "Start", "End", "Strand",
            "Length", "Source_File"
        ])
        writer.writerows(summary_rows)

    # pseudogenes summary
    pseudogene_rows = [row for row in summary_rows if row[0] == "pseudogene"]
    if pseudogene_rows:
        with open(output_dir / "pseudogenes" / "pseudogenes_summary.csv", "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([
                "Feature_Type", "Pseudogene_Name",
                "Start", "End", "Strand",
                "Length", "Source_File"
            ])
            writer.writerows(pseudogene_rows)

    copy_filtered_files(output_dir)
    combine_to_3_combined(output_dir)
    create_phylo_folder(output_dir)

    print(f"\n✅ Done! All results saved in: {output_dir}")

if __name__ == "__main__":
    main()

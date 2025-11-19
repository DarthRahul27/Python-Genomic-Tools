#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess

############################################################
# Parse: sample \t path/to/windows10kb/file
############################################################
def parse_sample_list(file_path):
    sample_files = {}
    with open(file_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            name, path = line.strip().split("\t")
            sample_files[name] = path
    return sample_files


############################################################
# Extract cumulative chromosome boundaries
############################################################
def get_chrom_boundaries(windows_file):
    chrom_order = []
    chrom_end = {}

    with open(windows_file) as fh:
        next(fh)
        for line in fh:
            p = line.strip().split("\t")
            chrom = p[0]
            run = int(p[3])

            if chrom not in chrom_end:
                chrom_order.append(chrom)
            chrom_end[chrom] = max(chrom_end.get(chrom, 0), run)

    boundaries = [chrom_end[c] for c in chrom_order]
    return boundaries, chrom_order


############################################################
# Generate base-R plot script with pastel colours + labels
############################################################
def generate_r_script(sample_files, features, y_min, y_max, output_file):

    pastel_cols = [
        "#A6CEE3", "#FDBF6F", "#B2DF8A",
        "#CAB2D6", "#FB9A99", "#FFFF99"
    ]

    with open(output_file, "w") as out:

        out.write("### AUTO-GENERATED R SCRIPT ###\n")
        out.write("pdf(file=\"%s.pdf\", width=14, height=10)\n\n" % output_file)

        # Colours and feature vector
        out.write("plot_cols <- c(%s)\n" %
                  ",".join([f'"{c}"' for c in pastel_cols]))
        out.write("features <- c(%s)\n\n" %
                  ",".join([f'"{f}"' for f in features]))

        # Load windows files
        out.write("df_list <- list()\n\n")
        for sample, path in sample_files.items():
            out.write(f"df <- read.table('{path}', sep='\\t', header=TRUE)\n")
            out.write(f"df$Sample <- '{sample}'\n")
            out.write("df_list[[length(df_list)+1]] <- df\n\n")
        out.write("all_df <- do.call(rbind, df_list)\n\n")

        # Chromosome boundaries
        first = next(iter(sample_files.values()))
        boundaries, chroms = get_chrom_boundaries(first)

        out.write("chrom_boundaries <- c(%s)\n" %
                  ",".join(map(str, boundaries)))
        out.write("chrom_names <- c(%s)\n\n" %
                  ",".join([f'"{c}"' for c in chroms]))

        # Layout
        out.write("samples <- unique(all_df$Sample)\n")
        out.write("par(mfrow=c(length(samples),1), mar=c(4,4,3,1))\n\n")

        #######################################################
        # PANEL LOOP
        #######################################################
        out.write("for (s in samples) {\n")
        out.write("  sub <- all_df[all_df$Sample == s,]\n")
        out.write("  Mb <- sub$Running_Position / 1e6\n")
        out.write("  ymin <- %f\n" % y_min)
        out.write("  ymax <- %f\n" % y_max)

        # Base empty plot setup
        out.write("  plot(Mb, sub[[features[1]]], type='n', ylim=c(ymin,ymax),\n")
        out.write("       xlab='Genome position (Mb)', ylab='', main=s)\n\n")

        #######################################################
        # CORRECT CHROMOSOME BOUNDARIES + LABEL POSITIONS
        #######################################################
        out.write("  chr_mb <- chrom_boundaries / 1e6\n")

        # Chromosome boundary lines
        out.write("  for (i in seq_along(chr_mb)) {\n")
        out.write("    abline(v = chr_mb[i], col='black', lwd=1.2, lty=3)\n")
        out.write("  }\n\n")

        # Correct chromosome midpoints
        out.write("  chr_start <- c(0, chr_mb[-length(chr_mb)])\n")
        out.write("  chr_end   <- chr_mb\n")
        out.write("  chr_mid   <- (chr_start + chr_end) / 2\n")

        # Label chromosomes at top of panel
        out.write("  text(chr_mid, ymax * 0.95, labels=chrom_names,\n")
        out.write("       cex=0.8, font=2)\n\n")

        #######################################################
        # FEATURE CURVES
        #######################################################
        out.write("  for (i in seq_along(features)) {\n")
        out.write("    y <- sub[[features[i]]]\n")
        out.write("    lines(Mb, y, col=plot_cols[i], lwd=1.3)\n")
        out.write("  }\n\n")

        # Legend
        out.write("  legend('topright', legend=features,\n")
        out.write("         col=plot_cols[1:length(features)], lwd=2, bty='n')\n")

        out.write("}\n\n")
        out.write("dev.off()\n")

    print(f"[INFO] R script written: {output_file}")


############################################################
# Run R script
############################################################
def run_r_script(r_script_path):
    print(f"[INFO] Running Rscript: {r_script_path}")
    try:
        subprocess.run(["Rscript", r_script_path], check=True)
        print(f"[INFO] PDF created: {r_script_path}.pdf")
    except Exception as e:
        print(f"[ERROR] Failed: {e}")


############################################################
# CLI
############################################################
def main():
    parser = argparse.ArgumentParser(
        description="Genome-wide multi-feature base-R plotting"
    )

    parser.add_argument("-n", "--name_tab_windows_file", required=True)
    parser.add_argument("-f", "--feature", required=True)
    parser.add_argument("-y", "--y_min_max", default="0,50")
    parser.add_argument("-o", "--output", default="windows-plot-across-genome.R")

    args = parser.parse_args()
    features = args.feature.split(",")
    y_min, y_max = map(float, args.y_min_max.split(","))
    sample_files = parse_sample_list(args.name_tab_windows_file)

    generate_r_script(sample_files, features, y_min, y_max, args.output)
    run_r_script(args.output)


if __name__ == "__main__":
    main()

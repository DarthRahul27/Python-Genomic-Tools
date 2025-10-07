# r.anand@exeter.ac.uk
# Module initialization for VCF/Fasta/GFF tools

try:
    # Works when running as package
    from . import vcf, fasta, gff
except ImportError:
    # Works when running as a standalone script
    import modules.vcf as vcf
    import modules.fasta as fasta
    import modules.gff as gff

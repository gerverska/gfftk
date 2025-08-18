from .gff import dict2gff3, gff2dict


def sanitize(args):
    """
    Sanitize GFF3 and FASTA files by converting them to a standardized format and back.

    This function parses GFF3 and FASTA files using `gff2dict` to convert them into a standardized dictionary format. It then converts the dictionary back to GFF3 format using `dict2gff3`, effectively sanitizing the input files.

    Parameters:
        args (argparse.Namespace): Command-line arguments containing paths to GFF3 and FASTA files, output file path, and optional debug and url_encode flags.

    Returns:
        None
    """
    # Handle combined GFF3+FASTA files
    fasta_input = args.fasta
    if not fasta_input:
        # Check if GFF3 file is a combined file
        from .gff import is_combined_gff_fasta

        if is_combined_gff_fasta(args.gff3):
            fasta_input = None  # Let gff2dict handle the combined file
        else:
            raise ValueError("FASTA file is required unless using a combined GFF3+FASTA file")

    Genes = gff2dict(args.gff3, fasta_input, debug=args.debug)
    dict2gff3(
        Genes,
        output=args.out,
        debug=args.debug,
        url_encode=getattr(args, "url_encode", False),
    )

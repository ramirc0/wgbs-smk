"""Clean SAM/BAM read names and add contig metadata to the `@SQ` header lines."""

import argparse
import sys

import pysam


def read_contig_file(filename: str) -> dict[str, dict[str, str]]:
    """Read contig metadata into a nested dict.

    Each line holds a contig name, then tab-separated `TAG:value` fields as in
    `@SQ` lines. Only `LN`, `M5`, `AS` and `SP` are kept.

    Parameters
    ----------
    filename : str
        Contig metadata file.

    Returns
    -------
    dict of str to dict of str to str
        Tags per contig name.
    """
    contigs = {}
    with open(filename) as f:
        for line in f:
            fields = line.strip().split("\t")
            if not fields:
                continue

            name = fields[0]
            tags = {}

            for field in fields[1:]:
                if field.startswith(("LN:", "M5:", "AS:", "SP:")):
                    tag_type = field[:2]
                    tag_value = field[3:]
                    tags[tag_type] = tag_value

            contigs[name] = tags

    return contigs


def clean_read_name(read_name: str) -> str:
    """Replace `@` and non-printable characters in a read name with `_`.

    Parameters
    ----------
    read_name : str
        Read name.

    Returns
    -------
    str
        Cleaned read name.
    """
    cleaned = []
    for char in read_name:
        if char == "@":
            cleaned.append("_")
        elif 33 <= ord(char) <= 126:
            cleaned.append(char)
        else:
            cleaned.append("_")

    return "".join(cleaned)


def update_header(
    header: pysam.AlignmentHeader, contig_metadata: dict[str, dict[str, str]] | None
) -> pysam.AlignmentHeader:
    """Add contig metadata tags to the `@SQ` header lines.

    Parameters
    ----------
    header : pysam.AlignmentHeader
        Input header.
    contig_metadata : dict of str to dict of str to str, optional
        Tags per contig name. The header is returned unchanged when empty.

    Returns
    -------
    pysam.AlignmentHeader
        Updated header.
    """
    if not contig_metadata:
        return header

    new_header = header.to_dict()

    if "SQ" in new_header:
        for sq_record in new_header["SQ"]:
            seq_name = sq_record.get("SN")
            if seq_name and seq_name in contig_metadata:
                for tag, value in contig_metadata[seq_name].items():
                    sq_record[tag] = value

    return pysam.AlignmentHeader.from_dict(new_header)


def build_parser():
    """Build the read_filter.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the read_filter.py flags.
    """
    parser = argparse.ArgumentParser(
        description="Filter SAM/BAM read names and update headers"
    )
    parser.add_argument("contig_file", nargs="?", help="Optional contig metadata file")
    parser.add_argument(
        "-i", "--input", default="-", help="Input SAM/BAM file (default: stdin)"
    )
    parser.add_argument(
        "-o", "--output", default="-", help="Output SAM/BAM file (default: stdout)"
    )
    return parser


def main(argv=None):
    """Copy the alignments with cleaned read names and the updated header.

    Writes SAM to stdout, or BAM to a named output file.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)

    contig_metadata = None
    if args.contig_file:
        contig_metadata = read_contig_file(args.contig_file)

    input_file = sys.stdin if args.input == "-" else args.input
    with pysam.AlignmentFile(input_file, "r") as infile:
        updated_header = update_header(infile.header, contig_metadata)

        output_file = sys.stdout if args.output == "-" else args.output
        output_mode = "w" if args.output == "-" else "wb"

        with pysam.AlignmentFile(
            output_file, output_mode, header=updated_header
        ) as outfile:
            for read in infile:
                read.query_name = clean_read_name(read.query_name)
                outfile.write(read)


if __name__ == "__main__":
    main()

"""
Unit tests for the gff module.
"""

import os
import tempfile

from gfftk.gff import (
    _detect_format,
    gff2dict,
    is_combined_gff_fasta,
    simplifyGO,
    split_combined_gff_fasta,
)

# pytest is not used in this file


class TestGFFParsing:
    """Tests for GFF parsing functions."""

    def test_detect_format(self):
        """Test the _detect_format function."""
        # Create a GFF3-like file for testing
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp:
            temp.write("##gff-version 3\n")
            temp.write("contig1\tprediction\tgene\t1\t1000\t.\t+\t.\tID=gene1;Name=test_gene\n")
            temp_name = temp.name

        try:
            # Test with a GFF3 file
            format_type = _detect_format(temp_name)
            # The function returns a tuple with a parser function and a format name
            assert isinstance(format_type, tuple)
            assert format_type[1] in ["gff3", "default"]
        finally:
            # Clean up
            os.unlink(temp_name)

        # Create a GTF-like file for testing
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp:
            temp.write(
                'contig1\tprediction\tgene\t1\t1000\t.\t+\t.\tgene_id "gene1"; transcript_id "mRNA1";\n'
            )
            temp_name = temp.name

        try:
            # Test with a GTF file
            format_type = _detect_format(temp_name)
            # The function returns a tuple with a parser function and a format name
            assert isinstance(format_type, tuple)
            assert format_type[1] in ["gtf", "default"]
        finally:
            # Clean up
            os.unlink(temp_name)

    def test_simplify_go(self):
        """Test the simplifyGO function."""
        # Test with a simple GO term
        go_terms = ["GO:0005634"]
        result = simplifyGO(go_terms)
        assert result == ["GO:0005634"]

        # Test with multiple GO terms
        go_terms = ["GO:0005634", "GO:0003677", "GO:0006355"]
        result = simplifyGO(go_terms)
        assert set(result) == set(["GO:0005634", "GO:0003677", "GO:0006355"])

        # Test with an empty list
        go_terms = []
        result = simplifyGO(go_terms)
        assert result == []


class TestCombinedGFFFormat:
    """Tests for combined GFF3+FASTA format functionality."""

    def test_is_combined_gff_fasta_true(self):
        """Test detection of combined GFF3+FASTA files."""
        # Create a combined file
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp:
            temp.write("##gff-version 3\n")
            temp.write("contig1\ttest\tgene\t1\t100\t.\t+\t.\tID=gene1\n")
            temp.write("##FASTA\n")
            temp.write(">contig1\n")
            temp.write("ATCGATCGATCG\n")
            temp_name = temp.name

        try:
            result = is_combined_gff_fasta(temp_name)
            assert result is True
        finally:
            os.unlink(temp_name)

    def test_is_combined_gff_fasta_false(self):
        """Test detection of regular GFF3 files."""
        # Create a regular GFF3 file
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp:
            temp.write("##gff-version 3\n")
            temp.write("contig1\ttest\tgene\t1\t100\t.\t+\t.\tID=gene1\n")
            temp_name = temp.name

        try:
            result = is_combined_gff_fasta(temp_name)
            assert result is False
        finally:
            os.unlink(temp_name)

    def test_split_combined_gff_fasta(self):
        """Test splitting combined GFF3+FASTA files."""
        # Create a combined file
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp:
            temp.write("##gff-version 3\n")
            temp.write("contig1\ttest\tgene\t1\t100\t.\t+\t.\tID=gene1\n")
            temp.write("##FASTA\n")
            temp.write(">contig1\n")
            temp.write("ATCGATCGATCG\n")
            temp_name = temp.name

        try:
            gff_content, fasta_content = split_combined_gff_fasta(temp_name)

            # Check GFF content
            gff_lines = gff_content.getvalue().split("\n")
            assert "##gff-version 3" in gff_lines
            assert "contig1\ttest\tgene\t1\t100\t.\t+\t.\tID=gene1" in gff_lines
            assert "##FASTA" not in gff_lines

            # Check FASTA content
            fasta_lines = fasta_content.getvalue().split("\n")
            assert ">contig1" in fasta_lines
            assert "ATCGATCGATCG" in fasta_lines

        finally:
            os.unlink(temp_name)


class TestNonStandardFeatures:
    """Tests for non-standard GFF3 features support."""

    def test_non_standard_features_parsing(self):
        """Test parsing of non-standard GFF3 features."""
        # Create a GFF3 file with non-standard features
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as gff_temp:
            gff_temp.write("##gff-version 3\n")
            gff_temp.write("contig1\ttest\tgene\t1\t1000\t.\t+\t.\tID=gene1\n")
            gff_temp.write("contig1\ttest\tmRNA\t1\t1000\t.\t+\t.\tID=mRNA1;Parent=gene1\n")
            gff_temp.write("contig1\ttest\texon\t1\t300\t.\t+\t.\tID=exon1;Parent=mRNA1\n")
            gff_temp.write("contig1\ttest\tintron\t301\t700\t.\t+\t.\tID=intron1;Parent=mRNA1\n")
            gff_temp.write("contig1\ttest\texon\t701\t1000\t.\t+\t.\tID=exon2;Parent=mRNA1\n")
            gff_temp.write(
                "contig1\ttest\tnoncoding_exon\t1\t100\t.\t+\t.\tID=nc_exon1;Parent=mRNA1\n"
            )
            gff_temp.write(
                "contig1\ttest\tfive_prime_UTR_intron\t101\t200\t.\t+\t.\tID=utr_intron1;Parent=mRNA1\n"
            )
            gff_temp.write(
                "contig1\ttest\tpseudogenic_exon\t1\t500\t.\t+\t.\tID=pseudo_exon1;Parent=gene1\n"
            )
            gff_name = gff_temp.name

        # Create a simple FASTA file
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as fasta_temp:
            fasta_temp.write(">contig1\n")
            fasta_temp.write("A" * 1000 + "\n")
            fasta_name = fasta_temp.name

        try:
            # Test that parsing doesn't fail with non-standard features
            result = gff2dict(gff_name, fasta_name, debug=False)

            # Should successfully parse the gene
            assert len(result) >= 1
            assert "gene1" in result

        finally:
            os.unlink(gff_name)
            os.unlink(fasta_name)

    def test_intron_processing(self):
        """Test processing of intron features to generate exon coordinates."""
        # Create a GFF3 file with introns
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as gff_temp:
            gff_temp.write("##gff-version 3\n")
            gff_temp.write("contig1\ttest\tgene\t100\t900\t.\t+\t.\tID=gene1;Name=gene1\n")
            gff_temp.write(
                "contig1\ttest\tmRNA\t100\t900\t.\t+\t.\tID=mRNA1;Parent=gene1;Name=mRNA1\n"
            )
            gff_temp.write("contig1\ttest\tintron\t301\t700\t.\t+\t.\tID=intron1;Parent=mRNA1\n")
            gff_name = gff_temp.name

        # Create a simple FASTA file
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as fasta_temp:
            fasta_temp.write(">contig1\n")
            fasta_temp.write("A" * 1000 + "\n")
            fasta_name = fasta_temp.name

        try:
            # Test that introns are processed to generate exon coordinates
            result = gff2dict(gff_name, fasta_name, debug=False)

            # Should successfully parse the gene
            assert len(result) >= 1
            assert "gene1" in result

            # Check that exon coordinates were generated correctly
            gene = result["gene1"]
            assert len(gene["mRNA"]) == 1
            assert len(gene["mRNA"][0]) == 2  # Should have 2 exons

            # Expected exons: (100, 300) and (701, 900)
            expected_exons = [(100, 300), (701, 900)]
            actual_exons = gene["mRNA"][0]
            assert actual_exons == expected_exons

        finally:
            os.unlink(gff_name)
            os.unlink(fasta_name)

    def test_multiple_introns_processing(self):
        """Test processing of multiple intron features."""
        # Create a GFF3 file with multiple introns
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as gff_temp:
            gff_temp.write("##gff-version 3\n")
            gff_temp.write("contig1\ttest\tgene\t100\t1000\t.\t+\t.\tID=gene1;Name=gene1\n")
            gff_temp.write(
                "contig1\ttest\tmRNA\t100\t1000\t.\t+\t.\tID=mRNA1;Parent=gene1;Name=mRNA1\n"
            )
            gff_temp.write("contig1\ttest\tintron\t201\t300\t.\t+\t.\tID=intron1;Parent=mRNA1\n")
            gff_temp.write("contig1\ttest\tintron\t501\t700\t.\t+\t.\tID=intron2;Parent=mRNA1\n")
            gff_temp.write("contig1\ttest\tintron\t801\t900\t.\t+\t.\tID=intron3;Parent=mRNA1\n")
            gff_name = gff_temp.name

        # Create a simple FASTA file
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as fasta_temp:
            fasta_temp.write(">contig1\n")
            fasta_temp.write("A" * 1000 + "\n")
            fasta_name = fasta_temp.name

        try:
            # Test that multiple introns are processed correctly
            result = gff2dict(gff_name, fasta_name, debug=False)

            # Should successfully parse the gene
            assert len(result) >= 1
            assert "gene1" in result

            # Check that exon coordinates were generated correctly
            gene = result["gene1"]
            assert len(gene["mRNA"]) == 1
            assert len(gene["mRNA"][0]) == 4  # Should have 4 exons

            # Expected exons: (100, 200), (301, 500), (701, 800), (901, 1000)
            expected_exons = [(100, 200), (301, 500), (701, 800), (901, 1000)]
            actual_exons = gene["mRNA"][0]
            assert actual_exons == expected_exons

        finally:
            os.unlink(gff_name)
            os.unlink(fasta_name)

    def test_non_standard_exon_types(self):
        """Test parsing of non-standard exon types."""
        # Create a GFF3 file with non-standard exon types
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as gff_temp:
            gff_temp.write("##gff-version 3\n")
            gff_temp.write("contig1\ttest\tgene\t100\t500\t.\t+\t.\tID=gene1;Name=gene1\n")
            gff_temp.write(
                "contig1\ttest\tmRNA\t100\t500\t.\t+\t.\tID=mRNA1;Parent=gene1;Name=mRNA1\n"
            )
            gff_temp.write(
                "contig1\ttest\tnoncoding_exon\t100\t200\t.\t+\t.\tID=nc_exon1;Parent=mRNA1\n"
            )
            gff_temp.write(
                "contig1\ttest\tpseudogenic_exon\t300\t400\t.\t+\t.\tID=pseudo_exon1;Parent=mRNA1\n"
            )
            gff_name = gff_temp.name

        # Create a simple FASTA file
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as fasta_temp:
            fasta_temp.write(">contig1\n")
            fasta_temp.write("A" * 1000 + "\n")
            fasta_name = fasta_temp.name

        try:
            # Test that non-standard exon types are parsed as regular exons
            result = gff2dict(gff_name, fasta_name, debug=False)

            # Should successfully parse the gene
            assert len(result) >= 1
            assert "gene1" in result

            # Check that exon coordinates were stored correctly
            gene = result["gene1"]
            assert len(gene["mRNA"]) == 1
            assert len(gene["mRNA"][0]) == 2  # Should have 2 exons

            # Expected exons: (100, 200), (300, 400)
            expected_exons = [(100, 200), (300, 400)]
            actual_exons = gene["mRNA"][0]
            assert actual_exons == expected_exons

        finally:
            os.unlink(gff_name)
            os.unlink(fasta_name)

    def test_sgd_gene_structure(self):
        """Test parsing of SGD-style gene structures with snoRNA_gene and snoRNA features."""
        # Create a GFF3 file with SGD gene structure
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as gff_temp:
            gff_temp.write("##gff-version 3\n")
            gff_temp.write(
                "chrI\tSGD\tsnoRNA_gene\t142367\t142468\t.\t+\t.\tID=YNCA0003W;Name=YNCA0003W\n"
            )
            gff_temp.write(
                "chrI\tSGD\tsnoRNA\t142367\t142468\t.\t+\t.\tID=YNCA0003W_snoRNA;Parent=YNCA0003W\n"
            )
            gff_temp.write(
                "chrI\tSGD\tnoncoding_exon\t142367\t142468\t.\t+\t.\tParent=YNCA0003W_snoRNA\n"
            )
            gff_name = gff_temp.name

        # Create a simple FASTA file
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as fasta_temp:
            fasta_temp.write(">chrI\n")
            fasta_temp.write("A" * 200000 + "\n")
            fasta_name = fasta_temp.name

        try:
            # Test that SGD gene structure is parsed correctly
            result = gff2dict(gff_name, fasta_name, debug=False)

            # Should successfully parse the gene
            assert len(result) >= 1
            assert "YNCA0003W" in result

            # Check that gene structure is correct
            gene = result["YNCA0003W"]
            assert gene["type"] == ["snoRNA"]
            assert gene["ids"] == ["YNCA0003W_snoRNA"]
            assert len(gene["mRNA"]) == 1
            assert len(gene["mRNA"][0]) == 1  # Should have 1 exon

            # Expected exon coordinates
            expected_exons = [(142367, 142468)]
            actual_exons = gene["mRNA"][0]
            assert actual_exons == expected_exons

        finally:
            os.unlink(gff_name)
            os.unlink(fasta_name)

"""Test code for the codons module."""

from codons import (
    split_codons, translate_codons, translate_dna
)


def test_split_codons() -> None:
    """Testing split_codons()."""
    assert split_codons("aaa") == ['aaa']
    assert split_codons("aaaccc") == ['aaa', 'ccc']
    assert split_codons("aaacccggg") == ['aaa', 'ccc', 'ggg']

    # Not a multiple of three
    assert split_codons("a") is None
    assert split_codons("ac") is None
    assert split_codons("acac") is None
    assert split_codons("acacg") is None

    # Empty string
    assert split_codons("") == []


def test_translate_codons() -> None:
    """Testing translate_codons()."""
    assert translate_codons(
        ["ATG", "ACC", "GAA", "CAG", "TAG"]
    ) == ['M', 'T', 'E', 'Q', '*']
    assert translate_codons(
        ["atg", "ACC", "gAa", "caG", "Tag"]
    ) == ['M', 'T', 'E', 'Q', '*']
    assert translate_codons(
        ["atg", "gAa", "Tag"]
    ) == ['M', 'E', '*']

    # Faulty codons
    assert translate_codons(["aaa", "cc", "ggg"]) is None
    assert translate_codons(["aaa", "foo", "ggg"]) is None
    assert translate_codons(["123", "ccc", "ggg"]) is None

    # Empty list
    assert translate_codons([]) == []


def test_translate_dna() -> None:
    """Testing translate_dna()."""
    assert translate_dna("ATGACCGAACAGTAG") == 'MTEQ*'
    assert translate_dna("atgaccgaacagtag") == 'MTEQ*'
    assert translate_dna("atggaatag") == 'ME*'

    # Not multiple of three
    assert translate_dna("aaaccggg") is None

    # Invalid codon
    assert translate_dna("aaaFOOggg") is None

    # Empty string
    assert translate_dna("") == ""

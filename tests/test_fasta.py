"""Proteome FASTA reading/validation rules (src/slparalog/fasta.py docstring)."""

import pytest

from slparalog.fasta import FastaError, canonical_order, clean_records, parse_fasta


def _read(text):
    return clean_records(parse_fasta(text.splitlines(True)))


def test_multiline_and_header_tokens():
    seqs, qc = _read(">g1 some description\nACDE\nFGHI\n>g2\nKLMN\n")
    assert seqs == {"g1": "ACDEFGHI", "g2": "KLMN"}
    assert qc.n_records == 2 and qc.n_kept == 2


def test_duplicate_ids_rejected():
    with pytest.raises(FastaError, match="duplicate"):
        _read(">g1\nACDE\n>g1 other\nACDF\n")


def test_terminal_stop_stripped_internal_stop_to_x():
    seqs, qc = _read(">g1\nACD*\n>g2\nAC*DE*\n")
    assert seqs == {"g1": "ACD", "g2": "ACXDE"}
    assert qc.n_terminal_stops_stripped == 2
    assert qc.n_internal_stops_to_x == 1


def test_nonstandard_letters_to_x_and_lowercase():
    seqs, qc = _read(">g1\nacdUBzjO\n")
    assert seqs["g1"] == "ACDXXXXX"
    assert qc.n_nonstandard_to_x == 5


def test_gaps_and_whitespace_removed():
    seqs, _ = _read(">g1\nAC-D. E\n")
    assert seqs["g1"] == "ACDE"


def test_empty_records_dropped_and_reported():
    seqs, qc = _read(">g1\n\n>g2\nACDE\n>g3\n*\n")
    assert list(seqs) == ["g2"]
    assert qc.dropped_empty == ["g1", "g3"]


def test_all_empty_is_error():
    with pytest.raises(FastaError):
        _read(">g1\n\n")


def test_illegal_characters_rejected():
    with pytest.raises(FastaError, match="illegal"):
        _read(">g1\nAC1DE\n")


def test_sequence_before_header_rejected():
    with pytest.raises(FastaError):
        _read("ACDE\n>g1\nAC\n")


def test_canonical_order_independent_of_input_order():
    a = {"x": "MKV", "y": "AAA", "z": "MKV"}
    b = {"z": "MKV", "x": "MKV", "y": "AAA"}
    assert canonical_order(a) == canonical_order(b) == ["y", "x", "z"]

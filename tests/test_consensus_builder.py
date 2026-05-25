## Test file 
import pytest

from genomekit.modules.consensus_builder import ConsensusBuilder

def test_cb_happy_path():
    sequences = ["ATGC", "ATGC", "ATGC"]

    cb = ConsensusBuilder(sequences)
    cb._validate()
    cb.build_msa()
    cb.build_profile()

    result = cb.consensus()

    assert cb.consensus() == "ATGC"

def test_cb_edge_case():
    sequences = ["ATGC"]

    cb = ConsensusBuilder(sequences)
    cb._validate()
    cb.build_msa()
    cb.build_profile()

    result = cb.consensus()

    assert result == "ATGC"

def test_error_case():
    try:
        cb = ConsensusBuilder([])
        cb._validate()
        assert False, "Expected ValueError for empty input"
    except ValueError:
        assert True

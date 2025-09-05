import json
from pathlib import Path

from ga4gh.va_spec.base.core import Statement

from biocommons.gks_conversion_tool.converter import convert_gks_to_hl7_v2


def test_convert_gks_to_hl7_v2():
    statement_json = json.load(Path.open("tests/statement.json"))
    statement = Statement(**statement_json)
    result = convert_gks_to_hl7_v2(statement)
    assert result is not None
    assert result["VARIANT_NAME"] == "EGFR T790M"
    assert result["CHROMOSOME"] == "NC_000007.13"
    assert result["GENE_STUDIED"] == "EGFR"
    assert result["DNA_CHANGE"] == "c.2369C>T"
    assert result["AMINO_ACID_CHANGE"] == "p.Thr790Met"
    assert result["GENOMIC_DNA_CHANGE"] == "g.55249071C>T"
    assert result["GENOMIC_REFERENCE_SEQUENCE_ID"] == "NC_000007.13"

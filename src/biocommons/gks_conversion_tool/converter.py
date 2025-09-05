"""Converter for GKS <-> HL7 v2"""

import logging
from typing import Any

from ga4gh.cat_vrs.models import Constraint
from ga4gh.va_spec.base.core import Statement
from ga4gh.vrs.models import Allele, Expression, SequenceLocation

_logger = logging.getLogger(__name__)

# TODO: make this a pydantic class to enforce required vs optional fields and types for the values
HL7V2 = {
    "VARIANT_NAME": "504",
    "DISCRETE_VARIANT": "505",
    # required - needed as 1, 2, 3, 4, X, Y, etc
    "CHROMOSOME": "510",
    "ALLELE_START": "511.1",
    "ALLELE_END": "511.2",
    "DNA_REGION": "513",
    # required
    "GENE_STUDIED": "514",
    "TRANSCRIPT_REFERENCE_SEQUENCE_ID": "516",
    # required
    "DNA_CHANGE": "518",
    # required
    "AMINO_ACID_CHANGE": "520",
    "MOLECULAR_CONSEQUENCE": "521",
    "PROTEIN_REFERENCE_SEQUENCE": "522",
    # required? - refseq for genomic change (ex: NC_000023.11)
    "GENOMIC_REFERENCE_SEQUENCE_ID": "524",
    # "AMPLIFICATION": "525", Not supporting in proof of concept - this is for CNVs and we are not handling those yet (only simply variants)
    "REFERENCE_ALLELE": "526",
    "OBSERVED_ALLELE": "527",
    # required
    "GENOMIC_DNA_CHANGE": "528",
    "CYTOGENETIC_LOCATION": "532",
    "PENETRANCE": "534",
    "GENETIC_VARIANT_SOURCE": "535",
    "ALLELE_LENGTH": "545",
    "STRUCTURAL_INNER_START": "546.1",
    "STRUCTURAL_INNER_END": "546.2",
    "STRUCTURAL_OUTER_START": "547.1",
    "STRUCTURAL_OUTER_END": "547.2",
    "COPY_NUMBER": "550",
    # "FUSED_GENES": "551", Not supported until Cat-VRS 2.0
    # hardcoded to "detected" for now - this is sample specific / unsure where to get from VA-Spec
    # required
    "VARIANT_ASSESSMENT": "552",
    # required
    "VARIANT_CLASSIFICATION": "553",
    "INTERPRETATION": "554",
    "MODE_OF_INHERITANCE": "560",
    # Experimental functional effect in va-spec
    "FUNCTIONAL_EFFECT": "561",
    "REPEAT_NUCLEOTIDES": "564",
    "REPEAT_NUMBER": "565",
    "AFFECTED_EXON_START": "572.1",
    "AFFECTED_EXON_END": "572.2",
    "AFFECTED_INTRON_START": "573.1",
    "AFFECTED_INTRON_END": "573.2",
    "INTERPRETATION_NOTE": "575",
}


def convert_gks_to_hl7_v2(statement: Statement) -> dict[str, Any]:
    """
    Convert a VA-Spec Statement to an HL7 v2-compatible dictionary of fields.

    Returns a dict keyed by HL7 field identifiers (see HL7V2 constants).
    Raises ValueError if required data are missing.
    """
    proposition = statement.proposition
    subject_variant = proposition.subjectVariant

    # 504 - Variant Name
    variant_name = subject_variant.name
    if not variant_name:
        _logger.warning("subjectVariant.name is missing or empty")
        # TODO: error here because I'm pretty sure this is required?
        variant_name = None

    # 505 - Discrete Genetic Variant (placeholder until models solidify)
    # TODO: need to wait for models for this or find out what expected format is

    constraints = subject_variant.constraints or []
    allele, location = None, None
    if constraints:
        allele, location = _find_genomic_allele_and_location(constraints)
    else:
        err = "subjectVariant.constraints is missing or empty"
        raise ValueError(err)

    # Get hgvs.g expression from the allele (e.g., 'NC_000007.13:g.140453136A>T')
    genomic_expression = _find_expression(allele, syntax="hgvs.g")
    hgvs_g = genomic_expression.value if genomic_expression else None
    # 524 - Genomic Reference Sequence ID
    chromosome_ref_seq, g_dot = _parse_hgvs_dot(hgvs_g)

    # 511 - Allele start/end
    # to get genomic allele, get the seqRef Id, give to seqrepo, convert to ncbi namespace, check prefix for NC or NG
    allele_start, allele_end = _get_location_interval(location)

    # 513 - DNA Region

    # 514 - Gene Studied
    gene_studied = proposition.geneContextQualifier.name

    # 516 - Transcript Reference Sequence ID

    # 518 - DNA Change
    coding_expression = _find_expression(allele, syntax="hgvs.c")
    hgvs_c = coding_expression.value if coding_expression else None
    c_dot = _parse_hgvs_dot(hgvs_c)[1]

    # 520 - Amino Acid Change
    protein_expression = _find_expression(allele, syntax="hgvs.p")
    hgvs_p = protein_expression.value if protein_expression else None
    p_dot = _parse_hgvs_dot(hgvs_p)[1]

    # 521 - Molecular Consequence - on hold until approved in va-spec

    # 522 - Protein Reference Sequence
    # get protein allele, get refGetAccession, use seqrepo to convert to ncbi namespace

    # 526 - Reference Allele

    # 527 - Observed Allele

    # 528 - Genomic DNA Change

    # 532 - Cytogenetic Location

    # 534 - Penetrance

    # 535 - Genetic Variant Source

    # 545 - Allele Length

    # 546 - Structural Inner Start/End

    # 547 - Structural Outer Start/End

    # 550 - Copy Number

    # 553 - Variant Classification
    variant_classification = statement.classification

    # 554 - Interpretation

    # 560 - Mode of Inheritance

    # 561 - Functional Effect

    # 564 - Repeat Nucleotides

    # 565 - Repeat Number

    # 572 - Affected Exon Start/End

    # 573 - Affected Intron Start/End

    # 575 - Interpretation Note

    result: dict[str, Any] = HL7V2.copy()  # start with all keys
    result["VARIANT_NAME"] = variant_name
    # TODO: this needs to be converted to shorthand
    result["CHROMOSOME"] = chromosome_ref_seq
    result["ALLELE_START"] = allele_start
    result["ALLELE_END"] = allele_end
    result["GENE_STUDIED"] = gene_studied
    result["DNA_CHANGE"] = c_dot
    result["AMINO_ACID_CHANGE"] = p_dot
    result["GENOMIC_DNA_CHANGE"] = g_dot
    result["GENOMIC_REFERENCE_SEQUENCE_ID"] = chromosome_ref_seq
    result["VARIANT_CLASSIFICATION"] = variant_classification

    return result


# --- Helpers: extract from VA objects -------------------------------------


def _find_genomic_allele_and_location(
    constraints: list[Constraint],
) -> tuple[Allele, SequenceLocation] | None:
    """
    From a list of constraints, return the first (allele, location)
    """
    for constraint in constraints:
        if constraint.root.type != "DefiningAlleleConstraint":
            continue
        allele = constraint.root.allele
        location = allele.location
        if location is None:
            continue
        return allele, location
    return None


def _find_expression(allele: Allele, syntax: str) -> Expression | None:
    """
    Find an expression with a given syntax (e.g., 'hgvs.g') from allele.expressions.
    Returns the first matching expression found.
    """
    expressions = allele.expressions or []

    for expr in expressions:
        s = expr.syntax
        if s == syntax:
            return expr
    # TODO: raise error?
    return None


def _get_location_interval(location: SequenceLocation) -> tuple[int, int]:
    """
    Extract (start, end) from a SequenceLocation.
    """
    start = location.start
    end = location.end
    return start, end


# --- Helpers: transformation / parsing ---------------------------------------


def _parse_hgvs_dot(hgvs_value: str) -> tuple[str, str]:
    """
    Parse an hgvs.(g,c,p) expression.

    Expected styles:
      - 'NC_000007.13:g.140453136A>T'

    Returns:
      (chromosome, dot) where chromosome is the left of ':', and dot includes g.,c.,or p. onwards.
    """
    chromosome, dot = hgvs_value.split(":", 1)

    return chromosome, dot

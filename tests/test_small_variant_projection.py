from backend.domain.small_variant_projection import (
    AssemblyBaseMap,
    ProjectedAllele,
    ReferenceBase,
    expected_transcript_reference,
    project_exact_indel,
)


TRANSCRIPT = "NM_007294.4"


def _assembly_map() -> AssemblyBaseMap:
    first = ReferenceBase("17", 101, "C", "C")
    last = ReferenceBase("17", 102, "G", "G")
    return AssemblyBaseMap(
        strand=1,
        by_transcript_position={(2, 0): first, (3, 0): last},
        by_genomic_position={100: "A", 101: "C", 102: "G"},
    )


def test_projection_and_reference_use_the_same_mapped_cds_interval():
    assembly_map = _assembly_map()

    assert expected_transcript_reference(
        TRANSCRIPT, "c.2_3del", assembly_map
    ) == "CG"
    assert project_exact_indel(
        "BRCA1", TRANSCRIPT, "c.2_3del", assembly_map
    ) == ProjectedAllele("17", 100, "ACG", "A")


def test_unmapped_interval_remains_unavailable():
    assembly_map = _assembly_map()

    assert expected_transcript_reference(TRANSCRIPT, "c.1_3del", assembly_map) is None
    assert project_exact_indel(
        "BRCA1", TRANSCRIPT, "c.1_3del", assembly_map
    ) is None

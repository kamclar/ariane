from scripts.verify_enigma_source_bundle import verify_bundle


def test_pinned_enigma_v1_2_source_bundle_is_intact():
    assert verify_bundle() == []

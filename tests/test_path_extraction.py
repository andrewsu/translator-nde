"""Unit tests for creative-mode path extraction.

Shaped after the Translator paper's use case 2 (celiprolol -> ADRB2 -> COL1A1 in
vEDS), with the mechanistic chain hidden behind an inferred `treats` edge's
`biolink:support_graphs` attribute -- which is how it actually arrives.
"""

from translator_nde.translator import (
    collect_edge_ids, dedupe_paths, extract_paths,
)

VEDS = "MONDO:0017314"

MESSAGE = {
    "knowledge_graph": {
        "nodes": {
            "CHEBI:3568": {"name": "celiprolol", "categories": ["biolink:SmallMolecule"]},
            "NCBIGene:154": {"name": "ADRB2", "categories": ["biolink:Gene"]},
            VEDS: {"name": "vascular Ehlers-Danlos", "categories": ["biolink:Disease"]},
        },
        "edges": {
            # The inferred summary edge the user first sees...
            "e_treats": {
                "subject": "CHEBI:3568", "object": VEDS, "predicate": "biolink:treats",
                "attributes": [
                    {"attribute_type_id": "biolink:support_graphs", "value": ["aux1"]}
                ],
            },
            # ...and the real chain, reachable only through the aux graph.
            "e_drug_gene": {
                "subject": "CHEBI:3568", "object": "NCBIGene:154",
                "predicate": "biolink:affects",
                "qualifiers": [
                    {"qualifier_type_id": "biolink:object_aspect_qualifier",
                     "qualifier_value": "activity"},
                    {"qualifier_type_id": "biolink:object_direction_qualifier",
                     "qualifier_value": "increased"},
                ],
                "sources": [{"resource_id": "infores:drugcentral",
                             "resource_role": "primary_knowledge_source"}],
            },
            "e_gene_disease": {
                "subject": "NCBIGene:154", "object": VEDS,
                "predicate": "biolink:gene_associated_with_condition",
            },
        },
    },
    "auxiliary_graphs": {
        # Deliberately reversed: the spec says this list is unordered.
        "aux1": {"edges": ["e_gene_disease", "e_drug_gene"]},
    },
    "results": [
        {
            "node_bindings": {"sn": [{"id": "CHEBI:3568"}], "on": [{"id": VEDS}]},
            "analyses": [
                {"edge_bindings": {"t_edge": [{"id": "e_treats"}]}, "score": 0.87}
            ],
        }
    ],
}


def test_follows_edge_attribute_support_graphs():
    """The chain hangs off an edge attribute, not analyses[].support_graphs."""
    ids = collect_edge_ids(MESSAGE, MESSAGE["results"][0])
    assert {"e_treats", "e_drug_gene", "e_gene_disease"} <= ids


def test_recovers_celiprolol_adrb2_path():
    paths = dedupe_paths(extract_paths(MESSAGE, VEDS))
    assert len(paths) == 1
    p = paths[0]
    assert (p.drug_name, p.gene_name) == ("celiprolol", "ADRB2")
    assert p.disease == VEDS
    assert p.score == 0.87


def test_captures_qualifiers_for_route_a_comparison():
    """Direction/aspect must survive -- they are what GXA is checked against."""
    p = dedupe_paths(extract_paths(MESSAGE, VEDS))[0]
    assert p.drug_gene.direction == "increased"
    assert p.drug_gene.aspect == "activity"
    assert p.drug_gene.sources == ["infores:drugcentral"]


def test_ignores_unordered_aux_edge_list():
    """Reversing the aux edge list must not change the extracted path."""
    import copy
    msg = copy.deepcopy(MESSAGE)
    msg["auxiliary_graphs"]["aux1"]["edges"] = ["e_drug_gene", "e_gene_disease"]
    assert [p.to_dict() for p in dedupe_paths(extract_paths(msg, VEDS))] == \
           [p.to_dict() for p in dedupe_paths(extract_paths(MESSAGE, VEDS))]


def test_rejects_non_gene_intermediate():
    """A drug -> disease -> disease chain must not be reported as drug->gene."""
    import copy
    msg = copy.deepcopy(MESSAGE)
    msg["knowledge_graph"]["nodes"]["NCBIGene:154"]["categories"] = ["biolink:Disease"]
    assert dedupe_paths(extract_paths(msg, VEDS)) == []


def test_biologic_is_not_a_gene_intermediate():
    """An antibody drug is biolink:Protein too -- it must not fill the gene slot.

    Regression for observed output `Etoricoxib -> Infliximab -> RA`.
    """
    import copy
    msg = copy.deepcopy(MESSAGE)
    msg["knowledge_graph"]["nodes"]["NCBIGene:154"]["categories"] = [
        "biolink:Protein", "biolink:Drug",
    ]
    assert dedupe_paths(extract_paths(msg, VEDS)) == []

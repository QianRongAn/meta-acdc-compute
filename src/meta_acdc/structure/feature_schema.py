"""Canonical graph layout. Appended ablations must preserve this prefix."""
AA_ORDER = 'ACDEFGHIKLMNPQRSTVWY'
AA_START = 4
AA_STOP = 24
PEPTIDE_FLAG_INDEX = 24
CONFIDENCE_INDEX = 3
NODE_DIM = 26
EDGE_DIM = 16
FEATURE_SCHEMA = "interface-heavy-centroid-v1-26x16"


def validate_node_features(rows):
    for row in rows:
        if len(row) not in (NODE_DIM, NODE_DIM + 2):
            raise ValueError("unknown node feature layout")
        if row[PEPTIDE_FLAG_INDEX] not in (0, 1):
            raise ValueError("peptide flag must be binary in raw graph")
        identity = row[AA_START:AA_STOP]
        if any(v not in (0, 1) for v in identity) or sum(identity) != 1:
            raise ValueError("canonical one-hot amino-acid identity required")

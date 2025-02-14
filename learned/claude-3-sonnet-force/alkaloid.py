"""
Classifies: CHEBI:22315 alkaloid
"""
ring_matches = [mol.GetSubstructMatches(ring_pattern) for ring_pattern in ring_patterns]
"""
Classifies: CHEBI:22315 alkaloid
"""
has_heterocyclic_ring = any(ring.IsHeterocyclic() for ring in ring_info.AtomRings())
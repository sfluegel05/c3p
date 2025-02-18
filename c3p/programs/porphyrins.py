"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: Porphyrins – natural pigments whose core structure is built from 
four pyrrole units fused by four methine groups forming a macrocyclic system.
"""

from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    A porphyrin is defined as a natural pigment whose core structure consists of
    four pyrrole-like (5-membered aromatic rings with one nitrogen) rings fused 
    by four methine (carbon) bridges.
    
    This implementation uses the molecule's ring information rather than a simple 
    SMARTS search. It counts all rings of size 5 that contain exactly one nitrogen,
    and then inspects the fused connectivity (shared bridging carbons) across those rings.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a porphyrin, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information using the SSSR (smallest set of smallest rings).
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # Each ring is a tuple of atom indices.
    
    # Identify “pyrrole rings”: 5-membered rings containing exactly one nitrogen.
    pyrrole_rings = []
    for ring in atom_rings:
        if len(ring) == 5:
            n_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if n_count == 1:
                # Store as a set of atom indices for ease of intersection analysis.
                pyrrole_rings.append(set(ring))
    
    if len(pyrrole_rings) < 4:
        return False, (f"Found only {len(pyrrole_rings)} pyrrole-like rings; "
                       "need at least 4 for a porphyrin core")
    
    # Determine bridging methine groups by finding carbons shared between pairs of pyrrole rings.
    bridging_carbons = set()
    n_rings = len(pyrrole_rings)
    for i in range(n_rings):
        for j in range(i+1, n_rings):
            common_atoms = pyrrole_rings[i].intersection(pyrrole_rings[j])
            for idx in common_atoms:
                # Consider only carbons as bridging atoms.
                if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6:
                    bridging_carbons.add(idx)
    
    if len(bridging_carbons) < 4:
        return False, (f"Only found {len(bridging_carbons)} bridging methine group(s); "
                       "expected at least 4 for a porphyrin core")
    
    return True, "Molecule has at least 4 pyrrole-like rings fused via at least 4 bridging methine groups, consistent with a porphyrin core"

# Example usage:
if __name__ == "__main__":
    # You can test one of the provided SMILES (for example, the last provided porphyrin)
    test_smiles = "c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)[nH]5)n4)[nH]3"
    result, reason = is_porphyrins(test_smiles)
    print("Result:", result)
    print("Reason:", reason)
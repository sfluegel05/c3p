"""
Classifies: CHEBI:23763 pyrroline
"""
"""
Classifies: CHEBI:35688 pyrroline
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline based on its SMILES string.
    A pyrroline is an organic heteromonocyclic compound with a dihydropyrrole structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrroline, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for dihydropyrrole ring pattern
    dihydropyrrole_pattern = Chem.MolFromSmarts("[NR1]1[CR1][CR1][CR1][CR1]1")
    if not mol.HasSubstructMatch(dihydropyrrole_pattern):
        return False, "No dihydropyrrole ring found"
    
    # Check for aromaticity (should be non-aromatic)
    if mol.GetAromaticRingInfo().AromaticRings:
        return False, "Pyrroline ring must be non-aromatic"
    
    # Check for other rings (should be monocyclic)
    rings = mol.GetRingInfo().AtomRings()
    if len(rings) > 1:
        return False, "Pyrroline should be a monocyclic compound"
    
    # Check for heteroatoms other than nitrogen
    allowed_atoms = [6, 7]  # C, N
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, "Pyrroline should only contain C and N atoms"
        elif atom.GetAtomicNum() == 7:
            # Check if nitrogen is in a 5-membered ring
            is_in_5ring = False
            for ring in rings:
                if atom.GetIdx() in ring and len(ring) == 5:
                    is_in_5ring = True
                    break
            if not is_in_5ring:
                return False, "Nitrogen should be part of a 5-membered ring"
    
    return True, "Contains a non-aromatic dihydropyrrole ring system"
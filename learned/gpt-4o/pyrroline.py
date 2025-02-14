"""
Classifies: CHEBI:23763 pyrroline
"""
from rdkit import Chem

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline based on its SMILES string.
    A pyrroline is defined as any organic heteromonocyclic compound with a structure based on a dihydropyrrole.

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
    
    # Extended patterns to cover more pyrroline variations
    pyrroline_patterns = [
        Chem.MolFromSmarts("C1CC=NC1"),   # 2-pyrroline (aromatic ring)
        Chem.MolFromSmarts("C1C=CCN1"),   # 3-pyrroline
        Chem.MolFromSmarts("N1C=CCC1"),   # Reversed 2-pyrroline
        Chem.MolFromSmarts("[nH]1cccn1"),  # Aromatic 2-pyrroline (e.g., pyrrolo-decoated)
        Chem.MolFromSmarts("[nX3]1[c,C][c,C][c,C][c,C]1"),  # Generic pyrroline with flexible atom type
    ]
    
    for pattern in pyrroline_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a dihydropyrrole (pyrroline) ring"
    
    return False, "No pyrroline structure found"
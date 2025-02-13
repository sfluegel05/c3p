"""
Classifies: CHEBI:23763 pyrroline
"""
from rdkit import Chem

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline based on its SMILES string.
    A pyrroline is defined as any organic heteromonocyclic compound with a
    structure based on a dihydropyrrole.

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
    
    # Patterns representing common pyrroline forms
    pyrroline_patterns = [
        Chem.MolFromSmarts("C1CC=NC1"),    # 2-Pyrroline
        Chem.MolFromSmarts("C1C=CCN1"),    # 3-Pyrroline
        Chem.MolFromSmarts("N1C=CCC1"),    # Another form
        Chem.MolFromSmarts("[nH]1cccn1"),  # Aromatic 2-pyrroline form
        Chem.MolFromSmarts("[nX3]1[c,C][c,C][c,C][c,C]1"), # Generic pyrroline
        Chem.MolFromSmarts("N1[CH2][CH2]C=C1"),  # 1-Pyrroline
        Chem.MolFromSmarts("C1C[NH]CC=1"),       # Substituted dihydropyrrole
        Chem.MolFromSmarts("C2=CC1=NCC=C1C=C2"), # Pyrroline in larger, fused systems
        Chem.MolFromSmarts("n1cccc1=O"),          # Pyrrolone which might indicate pyrroline
    ]
    
    for pattern in pyrroline_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a dihydropyrrole (pyrroline) ring"
    
    return False, "No pyrroline structure found"
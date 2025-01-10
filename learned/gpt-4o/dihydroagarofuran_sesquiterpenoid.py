"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
from rdkit import Chem

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    A dihydroagarofuran sesquiterpenoid is any sesquiterpenoid with a dihydroagarofuran skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroagarofuran sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Refine SMARTS pattern for a dihydroagarofuran skeleton; consider complexity, rings, and substituents.
    dihydroagarofuran_patterns = [
        Chem.MolFromSmarts("C1(C)OC2CCC3CC(C)(C(C2)O1)C3"),
        Chem.MolFromSmarts("C1(C)OC2CCC3CC(C(O2)C1)C3"),
        Chem.MolFromSmarts("C1C(C)(O)C2CCC3CC(C)(C(O1)C2)C3")
        # Additional patterns can be added based on diversity observed in dihydroagarofuran sesquiterpenoids
    ]

    # Check each pattern against the molecule
    for pattern in dihydroagarofuran_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains dihydroagarofuran skeleton"

    return False, "No dihydroagarofuran skeleton found in molecule"
"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
from rdkit import Chem

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    A dihydroagarofuran sesquiterpenoid is any sesquiterpenoid containing a dihydroagarofuran skeleton.

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

    # Refine SMARTS pattern for a dihydroagarofuran skeleton; consider complexity, rings, and substituents
    # Let's create a more general pattern that encompasses known structural elements of the skeleton.
    dihydroagarofuran_patterns = [
        Chem.MolFromSmarts("C1C[C@@H]2[C@H]3[C@@H](C1)C[C@@H](O2)CC3"),
        Chem.MolFromSmarts("C1C[C@H]2C[C@H]3[C@@H](C1)C[C@@H](O2)CC3"),
        # Additional patterns can be based on varying ring junction stereochemistry and common ester linkages
    ]

    # Check each pattern against the molecule
    for pattern in dihydroagarofuran_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains dihydroagarofuran skeleton"

    return False, "No dihydroagarofuran skeleton found in molecule"
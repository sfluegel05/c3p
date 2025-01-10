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

    # Update SMARTS pattern for the dihydroagarofuran skeleton
    # Reflects core bicyclic structure with additional stereochemistry and oxygens
    dihydroagarofuran_patterns = [
        Chem.MolFromSmarts("C1[C@H]2[C@@H](O1)C[C@H]3[C@@H]2CC3"),
        Chem.MolFromSmarts("C1[C@@H]2[C@H](O1)C[C@@H]3[C@H]2CC3"),
        # Include other patterns accounting for different observed stereochemistries
    ]

    # Check each pattern against the molecule
    for pattern in dihydroagarofuran_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains dihydroagarofuran skeleton"

    return False, "No dihydroagarofuran skeleton found in molecule"
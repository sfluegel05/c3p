"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
from rdkit import Chem

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Attempts to determine if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    Dihydroagarofuran sesquiterpenoids are complex structures with specific ring systems and functional groups.

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

    # Example structural pattern could include fused bicyclic ring system typical in these compounds
    # This is a simplification. Actual patterns would be more complex.
    # Note: Such patterns are placeholders and may not be precise for the class without expert derivation.
    dihydroagarofuran_core_pattern = Chem.MolFromSmarts("C1[C@@H]2[C@H](O)[C@H](O)C[C@@H]1C2")
    
    if not mol.HasSubstructMatch(dihydroagarofuran_core_pattern):
        return False, "Core structure of dihydroagarofuran not found"

    # Further checks could be added for the extensive ester groups or other characteristic features
    # This would require detailed knowledge of the variations in functional groups

    return True, "Molecule matches dihydroagarofuran sesquiterpenoid core structure"
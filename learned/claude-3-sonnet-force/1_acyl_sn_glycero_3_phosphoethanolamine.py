"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule belongs to the 1-acyl-sn-glycero-3-phosphoethanolamine class based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule belongs to the class, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the overall structural pattern
    pattern = Chem.MolFromSmarts("[C@@H](OP(=O)(OCCN)[O-])[CH2][OH].[CH2][OC](=O)[CX4,CX3]")
    if not mol.HasSubstructMatch(pattern):
        return False, "Molecule does not match the 1-acyl-sn-glycero-3-phosphoethanolamine structural pattern"

    return True, "Molecule belongs to the 1-acyl-sn-glycero-3-phosphoethanolamine class"
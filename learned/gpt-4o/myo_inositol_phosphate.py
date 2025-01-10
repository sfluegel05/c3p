"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
from rdkit import Chem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a myo-inositol phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for myo-inositol backbone (cyclohexane with 6 -OH)
    myo_inositol_pattern = Chem.MolFromSmarts("C1(CO)C(O)C(O)C(O)C(O)C1O")

    # Check for the myo-inositol backbone
    if not mol.HasSubstructMatch(myo_inositol_pattern):
        return False, "No myo-inositol backbone found"
    
    # Define SMARTS pattern for phosphate group
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(O)=O")

    # Check for presence of phosphate groups
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate groups found"
    
    return True, "Contains myo-inositol backbone with one or more phosphate groups"
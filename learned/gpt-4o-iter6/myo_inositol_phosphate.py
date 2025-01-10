"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
from rdkit import Chem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    A myo-inositol phosphate is an inositol phosphate in which the inositol component
    has myo-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        (bool, str): True if molecule is a myo-inositol phosphate with reason, else False with reason.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for myo-inositol core with flexible stereochemistry
    myo_inositol_pattern = Chem.MolFromSmarts("[C@@H]1(O)C(O)C(O)C(O)C(O)C1(O)")
    
    # Check for presence of myo-inositol core
    if not mol.HasSubstructMatch(myo_inositol_pattern):
        return False, "No myo-inositol core structure found"
    
    # Define SMARTS pattern for phosphate group
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(O)=O")
    
    # Ensure presence of phosphate groups
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate groups found"
    
    # Assume at least one phosphate must be directly attached to the inositol core
    inositol_with_phosphate = Chem.MolFromSmarts("[C@@H]1(OP(O)(O)=O)C(O)C(O)C(O)C(O)C1(O)")
    if not mol.HasSubstructMatch(inositol_with_phosphate):
        return False, "Phosphate groups not directly attached to myo-inositol core"

    return True, "Contains myo-inositol core structure with phosphate groups in characteristic locations"

# Example test
smiles_example = "O[C@@H]1[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H]1OP(O)(O)=O"
print(is_myo_inositol_phosphate(smiles_example))
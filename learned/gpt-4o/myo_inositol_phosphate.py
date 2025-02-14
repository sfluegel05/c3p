"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    A myo-inositol phosphate is an inositol structure with myo-configuration that has phosphate groups attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a myo-inositol phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a general cyclohexane with multiple hydroxyl groups which represents the myo-inositol core
    general_inositol_core = Chem.MolFromSmarts("C1(C)(O)[C](O)[C](O)[C](O)[C](O)[C]1O")
    if not mol.HasSubstructMatch(general_inositol_core):
        return False, "Molecule does not match the general cyclohexane hydroxylated structure"

    # Identify phosphate group presence
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O")
    phosphate_count = len(mol.GetSubstructMatches(phosphate_pattern))
    if phosphate_count == 0:
        return False, "No phosphate groups found"

    # Avoid structures with long carbon chains that might confound myo-inositol core
    if any(len(chain) > 6 for chain in Chem.rdmolops.GetMolFrags(mol, asMols=False)):
        return False, "Long carbon chain detected, possibly non-inositol component"
    
    # Check total number of hydroxyl and phosphate groups around the cyclohexane core
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4H](O)C")
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_count < 5 or (hydroxyl_count + phosphate_count) < 6:
        return False, "Insufficient hydroxyl/phosphate groups for a typical myo-inositol phosphate"
    
    return True, "Molecule contains inositol structure with appropriate number of hydroxyl and phosphate groups attached"
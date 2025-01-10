"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines if a molecule is a 1-phosphatidyl-1D-myo-inositol based on its SMILES string.
    This type of molecule has a phosphatidyl group attached to position 1 of 1D-myo-inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-phosphatidyl-1D-myo-inositol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")
    
    # SMARTS pattern for 1D-myo-inositol: a six-membered ring with hydroxyl groups in the right stereochemistry
    # Updated to ensure all hydroxyls are in correct stereochemistry
    inositol_pattern = Chem.MolFromSmarts("[C@H]1([C@@H]([C@H]([C@@H]([C@H]([C@H]1O)O)O)O)O)O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return (False, "1D-myo-inositol pattern not found or incorrect stereochemistry")

    # SMARTS pattern for a generic phosphatidyl group: focusing on glycerol with two esterified acyl chains and phosphate
    phosphatidyl_pattern = Chem.MolFromSmarts("OC(=O)[C@H](CO[P](=O)(O)O)OC(=O)C")
    if not mol.HasSubstructMatch(phosphatidyl_pattern):
        return (False, "Phosphatidyl group pattern not detected or incorrect")

    # Ensure the phosphatidyl group is bonded to position 1
    phosphorus_pattern = Chem.MolFromSmarts("P([O-])(=O)O")
    inositol_match_atoms = mol.GetSubstructMatch(inositol_pattern)
    phosphorus_match_atoms = mol.GetSubstructMatch(phosphorus_pattern)

    if not inositol_match_atoms or not phosphorus_match_atoms:
        return (False, "Phosphatidyl or inositol pattern not fully matched in the molecule")

    # Check connectivity
    inositol_atom_1 = inositol_match_atoms[0]
    phosphorus_atom = phosphorus_match_atoms[0]

    # Verify the connection via position 1 of the inositol moiety
    for neighbor in mol.GetAtomWithIdx(inositol_atom_1).GetNeighbors():
        if neighbor.GetIdx() == phosphorus_atom:
            return (True, "Valid 1-phosphatidyl-1D-myo-inositol detected")
    
    return (False, "Phosphatidyl group not attached to position 1 of inositol")
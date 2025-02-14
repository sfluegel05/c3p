"""
Classifies: CHEBI:28874 phosphatidylinositol
"""
"""
Classifies: Phosphatidylinositol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol based on its SMILES string.
    A phosphatidylinositol is a glycerophosphoinositol having one phosphatidyl group
    esterified to one of the hydroxy groups of inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylinositol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone: C-C-C with hydroxyl groups
    glycerol_pattern = Chem.MolFromSmarts("[CH2][CH][CH2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for two ester groups attached to glycerol carbons (sn-1 and sn-2)
    ester_pattern = Chem.MolFromSmarts("C(=O)O[CH]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, expected at least 2"

    # Check for phosphate group attached to glycerol (sn-3 position)
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O[CH]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group attached to glycerol"

    # Check for inositol ring attached to phosphate group
    inositol_pattern = Chem.MolFromSmarts("OC1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)C1")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring attached to phosphate group"

    # Optionally, check molecular weight range if needed
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for phosphatidylinositol"

    return True, "Molecule is a phosphatidylinositol"
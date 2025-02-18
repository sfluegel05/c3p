"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
"""
Classifies: phosphatidylinositol phosphate (PIP)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate based on its SMILES.
    These molecules have a glycerol backbone with two fatty acids, a phosphate group,
    and a myo-inositol ring with at least one phosphate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a PIP, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for glycerol backbone with two ester groups and one phosphate
    # Glycerol pattern: C-O-C(=O) (two esters) and C-O-P (phosphate)
    glycerol_pattern = Chem.MolFromSmarts("[CH2]-[CH](-[CH2])-[OX2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone detected"

    # Find phosphate connected to glycerol
    phosphate_glycerol = Chem.MolFromSmarts("[OX2]-P(=O)([OX2])-[OX2]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_glycerol)
    if not phosphate_matches:
        return False, "No phosphate group attached to glycerol"

    # Check for myo-inositol ring (cyclohexane with multiple hydroxyls and at least one phosphate)
    # Inositol pattern: six-membered ring with O's and at least one phosphate
    inositol_pattern = Chem.MolFromSmarts("[C]1-[C]-[C]-[C]-[C]-[C]-1")
    inositol_matches = mol.GetSubstructMatches(inositol_pattern)
    if not inositol_matches:
        return False, "No inositol ring found"

    # Check inositol has at least one phosphate group
    inositol_phosphate = Chem.MolFromSmarts("[C][O]P(=O)([O-])[OX2]")
    if not mol.HasSubstructMatch(inositol_phosphate):
        return False, "Inositol lacks phosphate groups"

    # Check two fatty acid esters on glycerol
    ester_pattern = Chem.MolFromSmarts("[CX4][OX2]C(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    return True, "Glycerol with two esters, phosphate-linked myo-inositol with phosphates"
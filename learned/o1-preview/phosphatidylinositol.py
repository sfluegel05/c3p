"""
Classifies: CHEBI:28874 phosphatidylinositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
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

    # Check for inositol ring
    # Inositol is a cyclohexane ring with hydroxyl groups on each carbon
    inositol_pattern = Chem.MolFromSmarts('C1(CO)C(CO)C(CO)C(CO)C1O')
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "Inositol ring not found"

    # Check for phosphate group connected to inositol
    # Phosphate group connected via oxygen to inositol
    phosphate_inositol_pattern = Chem.MolFromSmarts('C1(CO)C(CO)C(CO)C(CO)C1OP([O-])(=O)O')
    if not mol.HasSubstructMatch(phosphate_inositol_pattern):
        return False, "Phosphate group connected to inositol not found"

    # Check for glycerol backbone connected to phosphate
    glycerol_phosphate_pattern = Chem.MolFromSmarts('OCC(O)COP([O-])(=O)O')
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "Glycerol backbone connected via phosphate not found"

    # Check for ester linkages at sn-1 and sn-2 positions
    ester_bond_pattern = Chem.MolFromSmarts('C(=O)OCC(O)CO')
    ester_matches = mol.GetSubstructMatches(ester_bond_pattern)
    if len(ester_matches) < 2:
        return False, "Less than 2 ester bonds found at sn-1 and sn-2 positions"

    return True, "Molecule matches phosphatidylinositol structure"
"""
Classifies: CHEBI:28874 phosphatidylinositol
"""
"""
Classifies: CHEBI:18179 phosphatidylinositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol based on its SMILES string.
    A phosphatidylinositol consists of inositol connected via a phosphate group to a
    glycerol backbone which has two fatty acid chains attached via ester bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylinositol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for inositol (myo-inositol pattern with 6 hydroxyls)
    inositol_pattern = Chem.MolFromSmarts('[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@@H]1O')
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No myo-inositol core found"

    # Check for phosphate group connected to inositol
    phosphate_pattern = Chem.MolFromSmarts('[O]P(=O)([O])[O]')
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group present"

    # Verify phosphate is connected to glycerol backbone
    # Glycerol pattern: two ester groups and one phosphate
    glycerol_pattern = Chem.MolFromSmarts('[CH2](OC(=O)*)[CH2](OC(=O)*)COP(=O)(O)')
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone with two esters and phosphate not found"

    # Check for at least two ester groups (fatty acids)
    ester_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[CX3](=O)[OX2H0]')))
    if ester_count < 2:
        return False, f"Only {ester_count} ester groups found, need at least 2"

    return True, "Phosphatidyl group esterified to myo-inositol via phosphate linkage"
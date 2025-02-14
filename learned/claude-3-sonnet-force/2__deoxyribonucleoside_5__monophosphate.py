"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
"""
Classifies: CHEBI:18115 2'-deoxyribonucleoside 5'-monophosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_2__deoxyribonucleoside_5__monophosphate(smiles: str):
    """
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2'-deoxyribonucleoside 5'-monophosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for deoxyribose sugar
    deoxyribose_pattern = Chem.MolFromSmarts("[C@H]([C@@H]([C@H]([C@H](CO)O)O)O)[C@@H]1[CH]O1")
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No 2'-deoxyribose sugar found"

    # Check for phosphate group
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(=O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for nucleobase
    nucleobase_pattern = Chem.MolFromSmarts("a1a2a3a4a5a6")  # Any aromatic ring system
    if not mol.HasSubstructMatch(nucleobase_pattern):
        return False, "No nucleobase found"

    # Check for connection between deoxyribose and phosphate
    deoxyribose_phosphate_pattern = Chem.MolFromSmarts("[C@H]([C@@H]([C@H]([C@H](CO)O)O)O)[C@@H]1[CH]O1OP(O)(=O)O")
    if not mol.HasSubstructMatch(deoxyribose_phosphate_pattern):
        return False, "Deoxyribose and phosphate not connected"

    # Check for connection between deoxyribose and nucleobase
    deoxyribose_nucleobase_pattern = Chem.MolFromSmarts("[C@H]([C@@H]([C@H]([C@H](CO)O)O)O)[C@@H]1[CH]O1a1a2a3a4a5a6")
    if not mol.HasSubstructMatch(deoxyribose_nucleobase_pattern):
        return False, "Deoxyribose and nucleobase not connected"

    # Additional checks
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 15 or num_atoms > 50:
        return False, "Molecule size outside expected range"

    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 1 or num_rings > 3:
        return False, "Number of rings outside expected range"

    return True, "Contains 2'-deoxyribose sugar with 5'-phosphate group and a nucleobase"
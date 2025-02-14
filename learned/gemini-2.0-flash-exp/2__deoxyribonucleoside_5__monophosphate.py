"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2__deoxyribonucleoside_5__monophosphate(smiles: str):
    """
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate based on its SMILES string.
    A 2'-deoxyribonucleoside 5'-monophosphate consists of a 2-deoxyribose sugar, a nucleobase
    attached at the 1' position, and a phosphate group at the 5' position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2'-deoxyribonucleoside 5'-monophosphate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the 2'-deoxyribose pattern (5-membered ring with specific substitutions)
    # Using SMARTS to find the ring and check for substitution pattern of H in C2 and C1
    # Note: The 5' carbon is also included, and its connection to oxygen
    deoxyribose_pattern = Chem.MolFromSmarts("[C@H]1[CH2][C@H]([CH2])[C@@H](O)[C@@H]([OX2])O1")
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No 2'-deoxyribose sugar found"

    # Check for the 2' hydrogen
    deoxy_2_hydrogen = Chem.MolFromSmarts("[C@H]1[CH2][C@H]([H])[C@@H](O)[C@@H]([OX2])O1")
    if not mol.HasSubstructMatch(deoxy_2_hydrogen):
        return False, "No hydrogen found on the 2' carbon"

    # Define the phosphate group pattern, connected via oxygen to a carbon
    phosphate_pattern = Chem.MolFromSmarts("[OX2][PX4](=[OX1])([OX2])([OX2])")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"
    
    # Check that the phosphate is attached to C5 of the deoxyribose
    phosphate_connection_pattern = Chem.MolFromSmarts("[C@H]1[CH2][C@H]([CH2])[C@@H](O)[C@@H]([OX2][PX4](=[OX1])([OX2])([OX2]))O1")
    if not mol.HasSubstructMatch(phosphate_connection_pattern):
        return False, "Phosphate not at 5' position"
    
   # Check for nucleobase connected to the deoxyribose C1.
    nucleobase_connection_pattern = Chem.MolFromSmarts("[NX3;!H0][C@H]1[CH2][C@H]([CH2])[C@@H](O)[C@@H]([OX2])O1")
    if not mol.HasSubstructMatch(nucleobase_connection_pattern):
        return False, "No nucleobase connected to 1' position"

    return True, "Molecule contains a 2'-deoxyribose with a phosphate on the 5' position and a nucleobase at the 1' position"
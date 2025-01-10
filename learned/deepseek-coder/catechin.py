"""
Classifies: CHEBI:23053 catechin
"""
"""
Classifies: CHEBI:23053 catechin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    A catechin is a member of the class of hydroxyflavan that has a flavan-3-ol skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catechin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible flavan-3-ol core structure pattern
    flavan_3_ol_pattern = Chem.MolFromSmarts("[C@H]1[C@H](O)[C@@H](Oc2ccccc12)c1ccc(O)c(O)c1")
    if not mol.HasSubstructMatch(flavan_3_ol_pattern):
        return False, "No flavan-3-ol core structure found"

    # Check for hydroxyl groups at key positions (3, 5, 7, 3', 4')
    hydroxyl_positions = [3, 5, 7, 3, 4]  # Positions relative to the flavan-3-ol core
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 3:
        return False, f"Found {len(hydroxyl_matches)} hydroxyl groups, need at least 3"

    # Check for the presence of a benzopyran ring with possible substitutions
    benzopyran_pattern = Chem.MolFromSmarts("O1C2=C(C=CC=C2)C=C1")
    if not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "No benzopyran ring found"

    # Check for the presence of a benzene ring with possible substitutions
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(benzene_pattern):
        return False, "No benzene ring found"

    # Check molecular weight - catechins typically >200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for catechin"

    return True, "Contains flavan-3-ol core structure with hydroxyl groups at key positions"
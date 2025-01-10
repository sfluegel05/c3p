"""
Classifies: CHEBI:47916 flavonoid
"""
"""
Classifies: CHEBI:47916 flavonoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    A flavonoid is characterized by a 1-benzopyran core with an aryl substituent at position 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavonoid core pattern: 1-benzopyran with an aryl substituent at position 2
    flavonoid_core_pattern = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[c]2[c][c][c][c][c]2-[O]-[C]3=[C][C]=[C][C]=3")
    
    # Check if the molecule contains the flavonoid core pattern
    if not mol.HasSubstructMatch(flavonoid_core_pattern):
        return False, "No flavonoid core structure found (1-benzopyran with aryl substituent at position 2)"

    # Check for the presence of an aryl group at position 2 of the pyran ring
    aryl_pattern = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1")
    matches = mol.GetSubstructMatches(aryl_pattern)
    if len(matches) < 2:
        return False, "No aryl substituent found at position 2 of the pyran ring"

    # Check for the presence of hydroxyl groups, which are common in flavonoids
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 1:
        return False, "No hydroxyl groups found, which are common in flavonoids"

    # Check molecular weight - flavonoids typically have a molecular weight > 200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for a flavonoid"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 10:
        return False, "Too few carbons for a flavonoid"
    if o_count < 2:
        return False, "Too few oxygens for a flavonoid"

    return True, "Contains 1-benzopyran core with an aryl substituent at position 2, characteristic of flavonoids"
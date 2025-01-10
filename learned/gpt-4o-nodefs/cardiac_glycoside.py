"""
Classifies: CHEBI:83970 cardiac glycoside
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a cardiac glycoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Substructure pattern for a steroid backbone
    steroid_pattern = Chem.MolFromSmarts('C1CCC2CC3CC(C2C1)CCC4C3CCCC4')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid core found"

    # Check for glycosidic linkage to sugars
    glycoside_pattern = Chem.MolFromSmarts('OC[C@H]1O[C@@H]([C@H](O)[C@@H](O)[C@H]1O]')
    if not mol.HasSubstructMatch(glycoside_pattern):
        return False, "No glycosidic linkage found"

    # Check for lactone ring (e.g., butenolide)
    lactone_pattern = Chem.MolFromSmarts('O=C1CCOC1')
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring found"
    
    # Further verification can involve looking at hydroxyl groups on steroid, sugar counts, etc.
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and atom.GetTotalDegree() == 1)
    if hydroxyl_count < 3:
        return False, f"Insufficient hydroxyl groups in the molecule: {hydroxyl_count} found"
    
    return True, "The structure contains a steroid core with glycosidic linkage and a lactone ring characteristic of cardiac glycosides"

# Example usage
# result = is_cardiac_glycoside("SMILES_STRING_HERE")
# print(result)
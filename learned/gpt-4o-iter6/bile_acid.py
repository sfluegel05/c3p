"""
Classifies: CHEBI:3098 bile acid
"""
from rdkit import Chem

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    Bile acids are hydroxy-5beta-cholanic acids with specific functional groups and stereochemistry.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the 5beta-cholanic acid framework
    cholanic_pattern = Chem.MolFromSmarts('[C@]2([C@H]([C@H]1[C@@H](C(CCC3[C@H](CC(=O)O)C=C(C)C3)CC1)CC2=O)C(C)C)C')
    if not mol.HasSubstructMatch(cholanic_pattern):
        return False, "5beta-cholanic acid framework not found"

    # Look for the presence of hydroxyl groups at specific positions
    hydroxy_pattern = Chem.MolFromSmarts('[C@H](O)')
    hydroxyl_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxyl_matches) < 3:
        return False, "Insufficient number of hydroxy groups (less than three)"

    # Ensure configuration includes 5beta stereochemistry
    fivebeta_stereo = Chem.MolFromSmarts('[C@@H]1[C@H](CCC2[C@@H](CCC3(C4=CCC(=O)CC4)C3)[C@@H]12)C')
    if not mol.HasSubstructMatch(fivebeta_stereo):
        return False, "5beta stereochemistry not correctly found"

    # Check for a carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "Carboxylic acid group not found"

    return True, "Valid bile acid: contains steroid core with hydroxy groups, the 5beta stereochemistry, and a terminal carboxylic group"

# Example usage:
# result, reason = is_bile_acid("O[C@@H]1[C@]2([C@]...")
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_hydroxy_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy monocarboxylic acid, defined as having a hydroxy group beta to the carboxy group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3-hydroxy monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Find carboxylic acid groups
    carboxyl_pattern = Chem.MolFromSmarts('C(=O)[OH]')
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    
    if not carboxyl_matches:
        return False, "No carboxylic acid group found"

    # Pattern for 3-hydroxy monocarboxylic acid
    # Matches C(OH)-C-C(=O)OH with optional substitution
    pattern = Chem.MolFromSmarts('[#6;!$(C=O)](-[OX2H])[#6][#6](=O)[OH]')
    matches = mol.GetSubstructMatches(pattern)

    if matches:
        # Check that the hydroxy group is not part of another carboxylic acid
        for match in matches:
            hydroxy_c = mol.GetAtomWithIdx(match[0])
            if len([n for n in hydroxy_c.GetNeighbors() if n.GetAtomicNum()==8]) == 1:
                return True, "Found 3-hydroxy monocarboxylic acid pattern"

    # Try alternate pattern that allows for alpha substitution
    pattern2 = Chem.MolFromSmarts('[#6;!$(C=O)](-[OX2H])-[#6]-[#6](=O)[OH]')
    matches2 = mol.GetSubstructMatches(pattern2)

    if matches2:
        for match in matches2:
            hydroxy_c = mol.GetAtomWithIdx(match[0]) 
            if len([n for n in hydroxy_c.GetNeighbors() if n.GetAtomicNum()==8]) == 1:
                return True, "Found 3-hydroxy monocarboxylic acid pattern"

    return False, "No 3-hydroxy monocarboxylic acid pattern found"
# Pr=1.0
# Recall=0.9411764705882353
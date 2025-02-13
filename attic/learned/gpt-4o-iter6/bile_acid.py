"""
Classifies: CHEBI:3098 bile acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Steroid core pattern (cholanic acid structure) - C27 framework
    steroid_pattern = Chem.MolFromSmarts('[C@]12[C@]3([C@](C[C@@H]([C@]4([C@@]3(CC[C@]4([C@@H](CCC(O)=O)C)[H])[H])C)O)[H])')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Steroid core matching cholanic acid structure not found"
        
    # Check for presence of hydroxy groups
    hydroxy_pattern = Chem.MolFromSmarts('[OH]')
    num_hydroxy = len(mol.GetSubstructMatches(hydroxy_pattern))
    if num_hydroxy < 3:
        return False, f"Found {num_hydroxy} hydroxy groups, need at least 3"
    
    # Check for the 5beta configuration stereochemistry
    # Since this is stereospecific, consider exploring specific compounds for more patterns
    
    # 5beta-configuration example pattern (simplified, real match would use stereochemistry)
    fivebeta_pattern = Chem.MolFromSmarts('[C@@H]1CCC[C@@H]2[C@]1(CCC3[C@@]2(CCC3)C)C')
    if not mol.HasSubstructMatch(fivebeta_pattern):
        return False, "5beta configuration not found"
    
    # Presence of a terminal carboxylic acid on the side chain
    carboxylic_pattern = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "Terminal carboxylic acid group not found"

    return True, "Contains steroid core with hydroxy groups, appropriate stereochemistry, and terminal carboxylic acid"

# Example of some bile acids would be tested using:
# result, reason = is_bile_acid("O[C@@H]1[C@]2([C@]...")
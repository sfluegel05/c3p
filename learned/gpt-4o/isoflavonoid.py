"""
Classifies: CHEBI:50753 isoflavonoid
"""
from rdkit import Chem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is defined as a 1-benzopyran with an aryl substituent at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1-benzopyran structure: a benzene ring fused to a 6-membered ring with an oxygen
    benzopyran_pattern = Chem.MolFromSmarts("c1cc2ccocc2c1") # Adjusted pattern for benzopyran
    if benzopyran_pattern is None:
        return False, "Error in SMARTS pattern definition"

    if not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "No 1-benzopyran backbone found"
        
    # Pattern for aryl group at position 3 (adjacent to oxygen in the pyranyl ring)
    # Using additional rules to ensure the correct attachment
    aryl_group_pattern = Chem.MolFromSmarts("c1ccc(-c2c(oc3ccc(-c4ccccc4)cc3)c(=O)cc2)c1")  # More general aryl attachment
    if aryl_group_pattern is None:
        return False, "Error in SMARTS pattern definition for aryl group"

    if not mol.HasSubstructMatch(aryl_group_pattern):
        return False, "No aryl substituent at position 3 found"

    return True, "Contains 1-benzopyran with an aryl substituent at position 3"

# Testing a sample isoflavonoid SMILES
result, reason = is_isoflavonoid("COc1ccc(cc1)-c1coc2cc(O)cc(O)c2c1=O")
print(result, reason)
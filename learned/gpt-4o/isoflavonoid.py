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

    # Define a more flexible pattern for 1-benzopyran
    # Consider approximation of benzopyran as a fused bicyclic structure with oxygen in one ring
    benzopyran_pattern = Chem.MolFromSmarts("c1ccc2c(c1)occ2") # A representation of benzopyran system
    if not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "No 1-benzopyran backbone found"
        
    # Check for an aryl group (aromatic ring) directly connected to the benzopyran system
    # Detect any aromatic ring connected to the fused ring system
    aryl_group_pattern = Chem.MolFromSmarts("c1ccc(-c2c(oc3ccccc3)c(=O)cc2)c2ccccc2")  # Simplified to detect substitution
    if not mol.HasSubstructMatch(aryl_group_pattern):
        return False, "No aryl substituent at position 3 found"

    return True, "Contains 1-benzopyran with an aryl substituent at position 3"

# Testing a sample isoflavonoid SMILES
result, reason = is_isoflavonoid("COc1ccc(cc1)-c1coc2cc(O)cc(O)c2c1=O")
print(result, reason)
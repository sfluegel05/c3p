"""
Classifies: CHEBI:75874 para-terphenyl
"""
from rdkit import Chem

def is_para_terphenyl(smiles: str):
    """
    Determines if a molecule is a para-terphenyl based on its SMILES string.
    A para-terphenyl consists of a 1,4-diphenylbenzene skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a para-terphenyl, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for para-terphenyl: c1ccc(cc1)-c1ccc(cc1)-c1ccccc1
    para_terphenyl_pattern = Chem.MolFromSmarts("c1ccc(cc1)-c1ccc(cc1)-c1ccccc1")
    
    # Check if the molecule matches the para-terphenyl pattern
    if mol.HasSubstructMatch(para_terphenyl_pattern):
        return True, "Contains 1,4-diphenylbenzene skeleton"
    else:
        return False, "Does not contain 1,4-diphenylbenzene skeleton"

# Example usage
example_smiles = "c1ccc(cc1)-c1ccc(cc1)-c1ccccc1"  # 1,4-diphenylbenzene itself
result, reason = is_para_terphenyl(example_smiles)
print(result, reason)
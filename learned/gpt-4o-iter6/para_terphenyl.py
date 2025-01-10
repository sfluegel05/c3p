"""
Classifies: CHEBI:75874 para-terphenyl
"""
from rdkit import Chem

def is_para_terphenyl(smiles: str):
    """
    Determines if a molecule is a para-terphenyl derivative based on its SMILES string.
    A para-terphenyl consists of a 1,4-diphenylbenzene skeleton, potentially with substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule resembles a para-terphenyl structure, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more general SMARTS pattern for para-terphenyl
    # This pattern looks for three phenyl rings with 1,4-linkage, allowing for flexibility in connections
    para_terphenyl_pattern = Chem.MolFromSmarts("c1ccc(cc1)-c2c(-c3ccc(cc3))cccc2")
    
    # Check if the molecule matches the flexible para-terphenyl pattern
    if mol.HasSubstructMatch(para_terphenyl_pattern):
        return True, "Contains a flexible 1,4-diphenylbenzene skeleton"
    else:
        return False, "Does not contain a 1,4-diphenylbenzene skeleton"

# Example usage
example_smiles = "c1ccc(cc1)-c1ccc(cc1)-c1ccccc1"  # 1,4-diphenylbenzene itself
result, reason = is_para_terphenyl(example_smiles)
print(result, reason)
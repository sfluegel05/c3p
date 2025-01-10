"""
Classifies: CHEBI:39410 1,2,4-triazines
"""
"""
Classifies: CHEBI:38323 1,2,4-triazine
"""
from rdkit import Chem

def is_1_2_4_triazines(smiles: str):
    """
    Determines if a molecule is a 1,2,4-triazine based on its SMILES string.
    A 1,2,4-triazine has a six-membered ring with nitrogen atoms at positions 1, 2, and 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1,2,4-triazine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible 1,2,4-triazine core pattern
    # This pattern matches a six-membered ring with nitrogens at positions 1, 2, and 4
    # and allows for any type of bond and any substituents on the ring
    triazine_pattern = Chem.MolFromSmarts("[nH0]1[cH0][nH0][nH0][cH0][cH0]1")
    if not mol.HasSubstructMatch(triazine_pattern):
        return False, "No 1,2,4-triazine core found"

    # Verify that the matched pattern has nitrogens at positions 1, 2, and 4
    matches = mol.GetSubstructMatches(triazine_pattern)
    for match in matches:
        atoms = [mol.GetAtomWithIdx(idx) for idx in match]
        # Check if the first three atoms in the match are nitrogens
        if atoms[0].GetAtomicNum() == 7 and atoms[2].GetAtomicNum() == 7 and atoms[3].GetAtomicNum() == 7:
            return True, "Contains a 1,2,4-triazine core with nitrogens at positions 1, 2, and 4"

    return False, "No 1,2,4-triazine core with nitrogens at positions 1, 2, and 4 found"
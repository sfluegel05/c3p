"""
Classifies: CHEBI:3098 bile acid
"""
from rdkit import Chem

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    A bile acid generally has a steroid core with various hydroxyl or keto groups,
    and a carboxylic acid group in the side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a bile acid, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved steroid core pattern considering general connectivity of cyclopentane-perhydrophenanthrene
    steroid_core_pattern = Chem.MolFromSmarts(
        "[C@@]1(C)[C@H]2[C@H]3CC[C@@H](CC3)[C@]2(C)CC1"
    )
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No steroid core found"

    # Check for hydroxyl (-OH) and/or keto (C=O) groups
    # More inclusive pattern to broadly capture these functional groups
    hydroxyl_keto_pattern = Chem.MolFromSmarts("[OX2H1,C=O,C(=O)[O;H1,-]]")
    if len(mol.GetSubstructMatches(hydroxyl_keto_pattern)) == 0:
        return False, "No hydroxyl or keto groups found"

    # Look for a carboxylic acid group (COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O[H1,-]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    return True, "Contains a steroid core with substituents and a carboxylic acid group."
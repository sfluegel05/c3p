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

    # Refined steroid core pattern considering three six-membered rings fused with a five-membered ring
    steroid_core_pattern = Chem.MolFromSmarts("C1CCC2C3C(C1)CC4CCCC(C3)C24")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No steroid core found"

    # Check for hydroxyl and/or keto groups
    # Specific number is not too critical, general presence is sufficient
    hydroxyl_keto_pattern = Chem.MolFromSmarts("[OX2H1,C=O]")
    if not mol.HasSubstructMatch(hydroxyl_keto_pattern):
        return False, "No hydroxyl or keto groups found"

    # Look for a carboxylic acid group, critical for bile acids
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    return True, "Contains a steroid core with substituents and a carboxylic acid group."
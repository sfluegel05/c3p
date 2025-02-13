"""
Classifies: CHEBI:24128 furanocoumarin
"""
"""
Classifies: Furanocoumarin
"""
from rdkit import Chem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin consists of a furan ring fused with a coumarin structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Furan ring pattern
    furan_pattern = Chem.MolFromSmarts("c1occc1")

    # Coumarin framework (Benzopyrone), broader definition to include fusions
    coumarin_pattern = Chem.MolFromSmarts("c1ccc2c(c1)oc(=O)cc2")

    # General pattern for the furanocoumarin fusion
    furanocoumarin_pattern = Chem.MolFromSmarts("c1ccc2c(c1)oc(=O)cc2-c3occc3")

    # Check for coumarin structure
    if not mol.HasSubstructMatch(coumarin_pattern):
        return False, "No coumarin structure found"

    # Check for furan ring presence
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring found"

    # Check for the fusion of a furan ring with a coumarin structure
    if mol.HasSubstructMatch(furanocoumarin_pattern):
        return True, "Contains a furanocoumarin fusion structure"

    return False, "No recognized furanocoumarin fusion detected"
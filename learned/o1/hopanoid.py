"""
Classifies: CHEBI:51963 hopanoid
"""
"""
Classifies: CHEBI:27439 hopanoid
"""
from rdkit import Chem

def is_hopanoid(smiles: str):
    """
    Determines if a molecule is a hopanoid based on its SMILES string.
    A hopanoid is a triterpenoid based on a hopane skeleton (pentacyclic triterpene).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hopanoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the hopane skeleton without substituents and stereochemistry
    hopane_smiles = "C1CC2CCC3CCCC4CCCC5CCCCC5C4C3C2C1"  # Simplified hopane skeleton
    hopane_mol = Chem.MolFromSmiles(hopane_smiles)
    if hopane_mol is None:
        return None, "Invalid hopane skeleton SMILES"

    # Check for the hopane skeleton in the molecule
    if not mol.HasSubstructMatch(hopane_mol):
        return False, "No hopane skeleton found"

    return True, "Contains hopane skeleton characteristic of hopanoids"
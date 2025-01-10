"""
Classifies: CHEBI:29348 fatty amide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is characterized by an amide group attached to a long carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define pattern for amide group [CX3](=O)[NX3]
    amide_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3]")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage found"

    # Check for long carbon chain (at least 8 carbons in a continuous chain)
    carbon_chain_pattern = Chem.MolFromSmarts("[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No long carbon chain found"

    # If both patterns are found, it is likely to be a fatty amide
    return True, "Contains an amide linkage with a long carbon chain"
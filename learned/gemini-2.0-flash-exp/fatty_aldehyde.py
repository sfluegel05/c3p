"""
Classifies: CHEBI:35746 fatty aldehyde
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde has a terminal aldehyde group and a long carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty aldehyde, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one terminal aldehyde group (C=O, connected to one other C)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[#6]")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if not aldehyde_matches:
        return False, "No terminal aldehyde group found"

    # Check for a long alkyl carbon chain attached to the aldehyde carbon.
    # The chain must have a total of at least 6 carbons. 
    # The carbon directly attached to the aldehyde is included in this count.
    # This pattern looks for the aldehyde group's carbon, then 5 additional
    # carbons of type C (not aromatic) connected via single bonds.
    long_chain_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "Carbon chain too short or not a fatty acid derivative."


    # Check for number of carbons, must have at least 6
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 6:
        return False, "Too few carbons for a fatty aldehyde"

    return True, "Has a terminal aldehyde group and a long carbon chain"
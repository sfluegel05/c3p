"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    A long-chain fatty alcohol has a chain length ranging from C13 to C22 with a hydroxyl group attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a hydroxyl group (-OH) attached to an aliphatic carbon
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H]")  # Aliphatic carbon with OH
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)

    if not hydroxyl_matches:
        return False, "No aliphatic hydroxyl group (-OH) found on a carbon chain"

    # Utilize RDKit's function to identify the longest carbon chain
    longest_chain_length = rdMolDescriptors.CalcLongestAliphaticChain(mol)
    
    # Check if the longest aliphatic chain meets the criteria for a long chain
    if longest_chain_length >= 13 and longest_chain_length <= 22:
        return True, "Has a hydroxyl group and the carbon chain length is within range for long-chain fatty alcohol"
    else:
        return False, f"Carbon chain length is {longest_chain_length}, expected between 13 and 22"

__metadata__ = {   'chemical_class': {   'name': 'long-chain fatty alcohol',
                                         'definition': 'A fatty alcohol with a chain length ranging from C13 to C22.'},
                   'attempt': 0,
                   'success': True,
                   'best': True,
                   'message': None}
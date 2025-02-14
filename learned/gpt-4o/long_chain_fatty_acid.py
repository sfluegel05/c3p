"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acid based on its SMILES string.
    A long-chain fatty acid is defined as a fatty acid with a carbon chain length from C13 to C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify carboxylic acid group pattern and ensure it's not part of a ring system
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;R0]")  # Not in a ring
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No standalone carboxylic acid group found"

    # Compute total number of carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Large aliphatic chain pattern (open to limited branching)
    # Allows for the variety of bonding styles within chains: single '~', double '=', triple '#'
    large_chain_pattern = Chem.MolFromSmarts("[CH2,CH,CH3]~*" * 12)  # Matches long aliphatic chain segments
    num_of_matches = len(mol.GetSubstructMatches(large_chain_pattern))
    
    # Check if it matches the count requirements
    if not (13 <= carbon_count <= 22) or num_of_matches == 0:
        return False, f"Total carbon chain count {carbon_count} is not within the long-chain fatty acid range or insufficient chain match"

    return True, f"Valid long-chain fatty acid with {carbon_count} carbon atoms"
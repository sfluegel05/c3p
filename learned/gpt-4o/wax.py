"""
Classifies: CHEBI:73702 wax
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_wax(smiles: str):
    """
    Determines if a molecule is a wax based on its SMILES string.
    A wax is typically an ester formed between long-chain fatty acids and long-chain alcohols.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ester group (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group found"

    # Identify carbon chain lengths
    carbon_chains = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            chain_length = 0
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6:
                    chain_length += 1
            carbon_chains.append(chain_length)
    
    # Check if there are sufficiently long carbon chains
    long_chain_threshold = 12  # Minimum chain length for wax, could adjust as needed
    long_chains = [length for length in carbon_chains if length >= long_chain_threshold]
    if len(long_chains) < 2:
        return False, f"Insufficient long carbon chains, found {len(long_chains)} long chains"

    # Check for exactly one ester linkage
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Incorrect number of ester groups, found {len(ester_matches)}"

    return True, "Contains long-chain molecules with one ester linkage characteristic of waxes"
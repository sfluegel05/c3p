"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid (VLCFA) based on its SMILES string.
    VLCFAs are characterized by having a carbon chain longer than 20 carbons with a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a VLCFA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group -O=C(O) pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid found"
        
    # Count carbons in the longest continuous carbon chain
    longest_chain_length = max(len(chain) for chain in rdmolops.GetMolFrags(mol, asMols=False, sanitizeFrags=False) if isinstance(chain, tuple))
    if longest_chain_length <= 20:
        return False, f"Longest carbon chain is {longest_chain_length} carbons, VLCFAs require more than 20"

    return True, "Contains a very long carbon chain with carboxylic acid group"
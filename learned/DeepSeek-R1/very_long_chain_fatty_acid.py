"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid (chain length ≥ C22) based on its SMILES string.
    
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

    # Find carboxylic acid groups (-COOH)
    carb_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    matches = mol.GetSubstructMatches(carb_acid_pattern)
    
    # Must have exactly one carboxylic acid group
    if len(matches) != 1:
        return False, f"Found {len(matches)} carboxylic acid groups (need exactly 1)"
    
    # Get the carbonyl carbon and check if aliphatic
    carbonyl_carbon = matches[0][0]
    if mol.GetAtomWithIdx(carbonyl_carbon).GetIsAromatic():
        return False, "Carboxylic acid attached to aromatic ring"

    # Find alpha carbon (adjacent to carbonyl)
    alpha_carbons = [n for n in mol.GetAtomWithIdx(carbonyl_carbon).GetNeighbors() 
                    if n.GetAtomicNum() == 6 and n.GetIdx() != matches[0][1]]
    
    if not alpha_carbons:
        return False, "No alpha carbon found"
    
    # Calculate longest chain from alpha carbon using RDKit's built-in function
    chain_atoms = rdMolDescriptors._CalcCarbons(mol, alpha_carbons[0].GetIdx(), doubleBondTerminal=True)
    chain_length = len(chain_atoms) + 1  # Include carbonyl carbon
    
    if chain_length >= 22:
        return True, f"Main carbon chain length ({chain_length}) ≥ 22"
    else:
        return False, f"Main carbon chain length ({chain_length}) < 22"
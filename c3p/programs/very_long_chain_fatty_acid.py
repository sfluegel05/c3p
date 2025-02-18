"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    A very long-chain fatty acid is a fatty acid with a carbon chain greater than C22.
    Ultra-long-chain fatty acids are a subset with chain length > C27.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Identify the carbon chain by searching for a series of C atoms
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]([H])([H])([H,CH2])") # match CH2 or CH groups
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if not chain_matches:
        return False, "No carbon chain found"

    # Estimate carbon chain length by looking at the longest chain.
    # The number of carbons needs to be greater than 22.
    
    max_chain_length = 0
    for match in chain_matches:
        current_chain_length = 0
        visited = [False] * mol.GetNumAtoms()
        for atom_index in match:
            
            if not visited[atom_index]:
              
                current_chain_length +=1
                visited[atom_index] = True
                
                # Check for neighbor atoms that are also part of chain
                neighbors = mol.GetAtomWithIdx(atom_index).GetNeighbors()
                
                for neighbor in neighbors:
                    neighbor_idx = neighbor.GetIdx()
                    if mol.GetAtomWithIdx(neighbor_idx).GetSymbol() == 'C' and mol.GetAtomWithIdx(neighbor_idx).GetTotalNumHs()>= 1 and not visited[neighbor_idx] :
                        visited[neighbor_idx]= True
                        current_chain_length+=1
        max_chain_length = max(max_chain_length, current_chain_length)
    
    if max_chain_length > 27:
        return True, f"Ultra-long-chain fatty acid (C>{max_chain_length})"
    elif max_chain_length > 22:
         return True, f"Very long-chain fatty acid (C>{max_chain_length})"
    else:
       return False, f"Not a very long-chain fatty acid (C={max_chain_length})"
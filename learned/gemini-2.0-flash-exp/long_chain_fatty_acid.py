"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
"""
Classifies: long-chain fatty acid (C13 to C22)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acid (C13 to C22) based on its SMILES string.
    A long-chain fatty acid is defined as a molecule containing a carboxyl group and
    a linear carbon chain of length 13 to 22 carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-COOH, -COO-, or ester)
    acid_pattern = Chem.MolFromSmarts("C(=O)[OX1H0,OX2H1,OX2]") # check for -COOH, -COO-, and esters
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group or ester found"

    # Get the carbon of the carboxyl group
    carboxylic_carbon_matches = mol.GetSubstructMatches(acid_pattern)
    if not carboxylic_carbon_matches:
      return False, "No carboxylic acid carbon found"
    carboxylic_carbon_idx = carboxylic_carbon_matches[0][0] # get the first carbon idx
    carboxylic_atom = mol.GetAtomWithIdx(carboxylic_carbon_idx)


    # Find the alpha carbon (carbon directly attached to the carbonyl carbon)
    alpha_carbons = []
    for neighbor in carboxylic_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:
            alpha_carbons.append(neighbor)

    if not alpha_carbons:
        return False, "No alpha carbon found"


    def trace_linear_chain(start_carbon, visited=None, chain_count = 0, previous_carbon = None):
        if visited is None:
            visited = set()
        if start_carbon.GetIdx() in visited:
            return chain_count
        
        visited.add(start_carbon.GetIdx())
        chain_count += 1
        next_carbons = []
        for neighbor in start_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != carboxylic_carbon_idx and (previous_carbon is None or neighbor.GetIdx() != previous_carbon.GetIdx()):
                next_carbons.append(neighbor)
        
        if len(next_carbons) == 0:
          return chain_count
        elif len(next_carbons) > 1 :
          # Branching, so terminate path
          return chain_count
        else:
          return trace_linear_chain(next_carbons[0], visited, chain_count, start_carbon)


    max_chain_length = 0
    for alpha_carbon in alpha_carbons:
      chain_length = trace_linear_chain(alpha_carbon)
      max_chain_length = max(chain_length, max_chain_length)

    if max_chain_length < 13 or max_chain_length > 22:
      return False, f"Chain length is {max_chain_length}, not between 13 and 22 carbons"
    
    return True, f"Molecule is a long-chain fatty acid with a chain length of {max_chain_length} carbons."
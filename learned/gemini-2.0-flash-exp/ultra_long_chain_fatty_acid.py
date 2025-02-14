"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
from rdkit import Chem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid
    (chain length > C27) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ultra-long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "Molecule does not contain a carboxylic acid group."

    # Find the carboxyl carbon
    matches = mol.GetSubstructMatches(carboxyl_pattern)
    carboxyl_carbon_idx = matches[0][0] # Get the index of the carbon in C=O

    # find a carbon connected to the carboxyl carbon
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_carbon_idx)
    chain_start_idx = None
    for neighbor in carboxyl_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:
            chain_start_idx = neighbor.GetIdx()
            break

    if chain_start_idx is None:
        return False, "No carbon chain connected to carboxyl group"

    # Follow the carbon chain
    current_atom_idx = chain_start_idx
    previous_atom_idx = carboxyl_carbon_idx
    carbon_count = 1 # we have already started at one carbon
    
    while True:
      current_atom = mol.GetAtomWithIdx(current_atom_idx)
      next_atom_idx = None

      # Find neighbor that is carbon and is not the previous
      for neighbor in current_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != previous_atom_idx:
                if next_atom_idx is not None: # if there are more than 1 neighbor, break as we are not following a single chain
                    next_atom_idx = None
                    break
                next_atom_idx = neighbor.GetIdx()


      if next_atom_idx is None: # if we don't find any neighbor or multiple neighbors, we have reached the end of the chain
        break

      # Continue following
      carbon_count += 1
      previous_atom_idx = current_atom_idx
      current_atom_idx = next_atom_idx


    if carbon_count > 27:
        return True, f"Longest carbon chain has {carbon_count} carbons, which is > 27."
    else:
        return False, f"Longest carbon chain has {carbon_count} carbons, which is <= 27."
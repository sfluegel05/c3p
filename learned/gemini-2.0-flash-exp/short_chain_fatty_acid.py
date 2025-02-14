"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.
    A short-chain fatty acid is an aliphatic monocarboxylic acid with a chain length of less than C6
    (maximum of 5 carbons in the alkyl chain) and no non-hydrocarbon substituents other than
    the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for non-hydrocarbon atoms (excluding H, C, O)
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num != 1 and atomic_num != 6 and atomic_num != 8:
            return False, "Non-hydrocarbon substituent found"
            

    # Check for carboxylic acid group with attached alpha carbon
    carboxylic_acid_pattern = Chem.MolFromSmarts("[C]-[C](=O)O")
    matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not matches:
        return False, "No carboxylic acid group with attached alpha carbon found"


    # Get the carbon index connected to the carbonyl
    carbonyl_carbon_idx = matches[0][1] # index of the carbon in C(=O)O
    alpha_carbon_idx = matches[0][0] #index of carbon connected to the carbonyl
    
    
    def get_chain_length(mol, current_atom_idx, visited_atoms, chain_length):
        
        visited_atoms.add(current_atom_idx)
        max_len = chain_length
        
        neighbors = []
        for neighbor in mol.GetAtomWithIdx(current_atom_idx).GetNeighbors():
          neighbor_idx = neighbor.GetIdx()
          bond = mol.GetBondBetweenAtoms(current_atom_idx, neighbor_idx)
          if neighbor.GetAtomicNum() == 6 and neighbor_idx not in visited_atoms and (bond.GetBondType() == Chem.rdchem.BondType.SINGLE or bond.GetBondType() == Chem.rdchem.BondType.DOUBLE) :
              neighbors.append(neighbor_idx)

        if len(neighbors) == 0:
          return chain_length
        
        for neighbor_idx in neighbors:
            max_len = max(max_len, get_chain_length(mol,neighbor_idx, set(visited_atoms), chain_length+1))

        return max_len
            
    
    chain_length = get_chain_length(mol,alpha_carbon_idx, set(), 1)
      
    if chain_length > 5:
         return False, "More than 5 carbons in the alkyl chain"
   
    return True, "Short-chain fatty acid criteria met"
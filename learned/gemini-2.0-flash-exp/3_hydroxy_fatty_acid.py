"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    A 3-hydroxy fatty acid is a fatty acid with a hydroxyl group on the carbon
    in the beta- or 3-position relative to the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for a hydroxy group
    hydroxy_pattern = Chem.MolFromSmarts("[OH]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if not hydroxy_matches:
         return False, "No hydroxy group found"


    # Get the carbon of the carboxyl group and get the 3-position carbons
    for match in mol.GetSubstructMatches(carboxylic_acid_pattern):
        carboxyl_c_index = match[0]
        
        # Check for a carbon at the 3-position connected to the carboxyl group
        carboxyl_c_atom = mol.GetAtomWithIdx(carboxyl_c_index)
        
        if not carboxyl_c_atom:
            return False, "Carboxyl C not found"
        
        
        neighbor_atoms = carboxyl_c_atom.GetNeighbors()
        
        if not neighbor_atoms:
            return False, "Carboxyl carbon does not have neighbors"
        
        if len(neighbor_atoms) == 1:
            alpha_c_atom = neighbor_atoms[0] # the alpha carbon is the neighbor, need to check if it is a C
            if alpha_c_atom.GetSymbol() != 'C':
                 return False, "Alpha neighbor must be a carbon"

        elif len(neighbor_atoms) > 1:
            alpha_c_atom = None
            for n_atom in neighbor_atoms:
               if n_atom.GetSymbol() == 'C' and n_atom.GetIdx() != carboxyl_c_index:
                  alpha_c_atom = n_atom
                  break
            if not alpha_c_atom:
                 return False, "No alpha carbon found"
        else:
            return False, "Cannot find alpha carbon"
            
        if not alpha_c_atom:
            return False, "Cannot determine alpha carbon"
    
        beta_c_neighbors = alpha_c_atom.GetNeighbors()
        
        if not beta_c_neighbors:
             return False, "Cannot find beta neighbors"

        if len(beta_c_neighbors) == 1: # end of chain
            return False, "Beta carbon at end of chain, does not fit"

        beta_c_atom = None
        for n_atom in beta_c_neighbors:
             if n_atom.GetSymbol() == 'C' and n_atom.GetIdx() != carboxyl_c_index and n_atom.GetIdx() != alpha_c_atom.GetIdx():
                 beta_c_atom = n_atom
                 break
        if not beta_c_atom:
            return False, "Cannot find beta carbon"
        
        gamma_c_neighbors = beta_c_atom.GetNeighbors()
        if not gamma_c_neighbors:
            return False, "Cannot find gamma neighbors"
    
        gamma_c_atom = None
        for n_atom in gamma_c_neighbors:
            if n_atom.GetSymbol() == 'C' and n_atom.GetIdx() != carboxyl_c_index and n_atom.GetIdx() != alpha_c_atom.GetIdx() and n_atom.GetIdx() != beta_c_atom.GetIdx():
                gamma_c_atom = n_atom
                break

        if not gamma_c_atom:
             return False, "Cannot find gamma carbon"

        for hydroxy_match in hydroxy_matches:
            hydroxy_atom = mol.GetAtomWithIdx(hydroxy_match[0])
            
            if hydroxy_atom.GetNeighbors():
              neighbor_of_hydroxy = hydroxy_atom.GetNeighbors()[0]
              
              if neighbor_of_hydroxy.GetIdx() == gamma_c_atom.GetIdx():
                break # found a gamma hydroxy group, we can continue the check.
        else:
           return False, "No hydroxy group found at gamma position"

    # Check fatty acid chain length (simplified to check for at least 4 carbons)
    chain_pattern = Chem.MolFromSmarts("C-C-C-C")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Too short hydrocarbon chain"

    return True, "3-Hydroxy fatty acid identified"
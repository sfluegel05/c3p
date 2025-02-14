"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an omega-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxyl group (C(=O)O or C(O)=O)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        carboxyl_pattern = Chem.MolFromSmarts("C(O)=O")
        if not mol.HasSubstructMatch(carboxyl_pattern):
            return False, "No carboxyl group found"

    # Check for at least one hydroxyl group (OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
         return False, "No hydroxyl group found"

    # Find position of the carbonyl carbon in the COOH
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    carbonyl_atom_idx = carboxyl_matches[0][0] if carboxyl_matches else None
    if carbonyl_atom_idx is None:
       return False, "Could not determine carbonyl position"

    # find all hydroxyl groups
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)

    # Get atoms
    atoms = mol.GetAtoms()
    n_atoms = len(atoms)
    if n_atoms < 3:
        return False, "Too few atoms to be an omega-hydroxy fatty acid"
    
    
    # Count total number of carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 3:
        return False, "Too few carbon atoms"
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
         return False, "Molecular weight too low"

    valid_omega = False
    reason = "No omega-hydroxyl group found"


    for match in hydroxyl_matches:
        hydroxyl_atom_idx = match[0]

        if hydroxyl_atom_idx == carbonyl_atom_idx:
           continue #same atom

        #Get the carbon chain
        
        #Find the atom connected to the carbonyl carbon, this will be the first C in chain
        carbonyl_atom = atoms[carbonyl_atom_idx]
        carbonyl_neighbors = [x.GetIdx() for x in carbonyl_atom.GetNeighbors() if x.GetAtomicNum() == 6]
        if len(carbonyl_neighbors) != 1:
           continue #This is not a linear chain if there is not one carbon atom adjacent to the carboxyl carbon
        first_chain_carbon_idx = carbonyl_neighbors[0]

        #Now the first atom in the chain is defined. From here, follow path to hydroxyl group
        
        current_atom_idx = first_chain_carbon_idx
        prev_atom_idx = carbonyl_atom_idx
        
        chain_atoms_list = []
        chain_atoms_list.append(carbonyl_atom_idx)
        chain_atoms_list.append(current_atom_idx)
        is_linear = True

        while current_atom_idx != hydroxyl_atom_idx:
              current_atom = atoms[current_atom_idx]
              
              next_carbons = [x.GetIdx() for x in current_atom.GetNeighbors() if x.GetAtomicNum() == 6 and x.GetIdx() != prev_atom_idx]
              if len(next_carbons) != 1:
                   is_linear = False
                   break
              prev_atom_idx = current_atom_idx
              current_atom_idx = next_carbons[0]

              chain_atoms_list.append(current_atom_idx)

        if is_linear:
            last_chain_atom = atoms[chain_atoms_list[-1]]
            
            #Make sure last chain atom has a hydroxyl group attached to it.
            
            found_hydroxy = False
            for neighbor in last_chain_atom.GetNeighbors():
                 if neighbor.GetAtomicNum() == 8:
                      found_hydroxy = True
                      break
            
            if found_hydroxy:

                  
                #Ensure that there are no other branches on chain.

                is_branching = False
                for i, idx in enumerate(chain_atoms_list):
                      
                      atom = atoms[idx]
                      
                      if i != 0 and i != len(chain_atoms_list)-1:
                        
                            carbon_neighbors = [x.GetIdx() for x in atom.GetNeighbors() if x.GetAtomicNum() == 6 and x.GetIdx() != chain_atoms_list[i-1] and x.GetIdx() != chain_atoms_list[i+1]]
                            if len(carbon_neighbors) > 0:
                                is_branching= True
                                break
                if not is_branching:
                  valid_omega = True
                  reason = "Omega hydroxyl group found on a linear chain"
                  break


    return valid_omega, reason
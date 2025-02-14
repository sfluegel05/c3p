"""
Classifies: CHEBI:87658 decanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester has a decanoyl group (10 carbons) attached via an ester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the ester group
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")

    # Find all ester groups
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    if not ester_matches:
        return False, "No ester bond found"
    
    # Check for the decanoyl group (10 carbon chain attached to C=O) for each ester
    for match in ester_matches:
        carbonyl_carbon = match[0] # C of C=O

        # Get the neighbor atoms of the carbonyl carbon
        neighbors = mol.GetAtomWithIdx(carbonyl_carbon).GetNeighbors()
        
        # Find the neighbor that is NOT the ester oxygen
        chain_start_atom = None
        for neighbor in neighbors:
            if neighbor.GetIdx() != match[1]:
                chain_start_atom = neighbor.GetIdx()
                break
        
        if chain_start_atom is None:
            continue # this is not the ester bond we need

        # get the atom count of the chain
        
        chain_atoms = []
        
        def find_chain(atom_index, current_chain,visited):
            visited.add(atom_index)
            current_chain.append(atom_index)
            
            atom = mol.GetAtomWithIdx(atom_index)
            neighbors = atom.GetNeighbors()
            for neighbor in neighbors:
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited and neighbor.GetAtomicNum() == 6:
                     find_chain(neighbor_idx, current_chain, visited)
            return current_chain
            
        chain_atoms = find_chain(chain_start_atom, [], set())

        if len(chain_atoms) == 10: # count the carbons in the chain
              
              
            # Now check that it is linear (no branches); if one of the atoms has more than 2 neighbors, then it is branched
            linear = True
            for idx in chain_atoms:
                if len(mol.GetAtomWithIdx(idx).GetNeighbors()) > 2 and idx != chain_start_atom and idx != chain_atoms[-1]:
                    linear = False
                    break
            if linear:
                  
                  
                # Verify it terminates with a -CH3 group
                end_atom_index = chain_atoms[-1] # the last atom of the chain
                
                end_atom = mol.GetAtomWithIdx(end_atom_index)
                
                if len(end_atom.GetNeighbors()) == 1 and end_atom.GetTotalNumHs() == 3: # the last atom should have 3 hydrogens
                       return True, "Contains a decanoyl group connected via an ester bond"
    
    
    return False, "No decanoate ester core found"
"""
Classifies: CHEBI:24026 fatty alcohol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.
    A fatty alcohol is defined as an aliphatic alcohol consisting of a chain of 3 to greater than 27 carbon atoms.
    Fatty alcohols may be saturated or unsaturated and may be branched or unbranched.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty alcohol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one alcohol group
    alcohol_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "No alcohol group found"
    
    # Initialize
    max_chain_length = 0

    # Find all alcohol groups
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)

    for match in alcohol_matches:
         # Get the carbon adjacent to the alcohol
        alcohol_atom = mol.GetAtomWithIdx(match[0])
        for neighbor in alcohol_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
              current_chain_length = 0
              visited_atoms = set()
              
              #Recursive function
              def traverse_chain(current_atom, chain_length, visited_atoms):
                nonlocal max_chain_length
                
                if current_atom.GetIdx() in visited_atoms:
                  return

                visited_atoms.add(current_atom.GetIdx())
                
                is_ring = False
                for bond in current_atom.GetBonds():
                    if bond.IsInRing():
                        is_ring = True
                        break
                if is_ring:
                    return
                
                chain_length +=1 
                max_chain_length = max(max_chain_length,chain_length)


                # Iterate through neighbor atoms 
                for neighbor in current_atom.GetNeighbors():
                  if neighbor.GetAtomicNum() == 6:
                    traverse_chain(neighbor,chain_length, visited_atoms)
              
              # start recursive chain
              traverse_chain(neighbor, 0, visited_atoms)
    
    # Check if a chain of appropriate length is found
    if max_chain_length < 3:
          return False, f"Too few carbons in a chain ({max_chain_length}), must be at least 3"

    #Check that there is at least one carbon in a chain that is not in a ring and that the number of carbon atoms is greater than 2
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    num_carbons = len(carbon_atoms)

    # Check if number of carbon is within the range or over
    if num_carbons < 3:
        return False, f"Too few carbon atoms ({num_carbons}), must be at least 3"
    
    return True, "Meets criteria for a fatty alcohol (3 to >27 C, at least one OH)"
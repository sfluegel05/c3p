"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: CHEBI:90255 epoxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.rdchem import Mol

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    An epoxy fatty acid contains a carboxylic acid group and an epoxide ring (3-membered cyclic ether),
    with the carboxylic acid being part of a carbon chain of at least 12 carbons.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an epoxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for carboxylic acid group (COOH or COO-)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1,O-]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxylic acid group"
    
    # Check for epoxide ring (3-membered cyclic ether)
    epoxide_pattern = Chem.MolFromSmarts("[O]1CC1")
    if not mol.HasSubstructMatch(epoxide_pattern):
        return False, "No epoxide ring found"
    
    # Check chain length for each carboxylic acid group
    max_chain_length = 0
    ring_info = mol.GetRingInfo()
    
    for match in carboxyl_matches:
        acid_carbon = match[0]  # Carbon in CX3(=O)
        visited = set()
        stack = [(acid_carbon, 1)]  # (atom index, current chain length)
        current_max = 0
        
        while stack:
            atom_idx, length = stack.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            
            if length > current_max:
                current_max = length
            
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx in visited:
                    continue
                
                # Check if neighbor is part of an epoxide ring
                if neighbor.GetAtomicNum() == 8:  # Oxygen
                    in_epoxide = False
                    for ring in ring_info.AtomRings():
                        if len(ring) == 3 and neighbor_idx in ring:
                            carbons = [a for a in ring if mol.GetAtomWithIdx(a).GetAtomicNum() == 6]
                            if len(carbons) == 2:
                                # Find the other carbon in the ring not already visited
                                other_carbon = next((a for a in ring if a != atom_idx and mol.GetAtomWithIdx(a).GetAtomicNum() == 6), None)
                                if other_carbon is not None and other_carbon not in visited:
                                    stack.append((other_carbon, length + 1))
                                    in_epoxide = True
                                    break
                    if in_epoxide:
                        continue  # Skip adding the oxygen to the stack, already handled
                
                # Traverse through non-oxygen atoms or non-epoxide oxygens
                if neighbor.GetAtomicNum() == 6:  # Carbon
                    stack.append((neighbor_idx, length + 1))
        
        if current_max > max_chain_length:
            max_chain_length = current_max
    
    if max_chain_length < 12:
        return False, f"Main chain too short ({max_chain_length} carbons), need ≥12"
    
    return True, "Contains carboxylic acid and 3-membered epoxide ring with sufficient chain length"
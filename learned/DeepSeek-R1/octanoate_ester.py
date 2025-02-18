"""
Classifies: CHEBI:87657 octanoate ester
"""
"""
Classifies: CHEBI:52077 octanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester has at least one ester group where the acid component is octanoic acid (CH3(CH2)6COO-).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octanoate ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Find all ester groups (O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"

    # Check each ester group for octanoate chain
    for match in ester_matches:
        oxygen_idx = match[0]
        carbonyl_carbon = match[1]

        # Get the adjacent atom to carbonyl (should be the acid's alkyl chain)
        acid_chain_atom = None
        for neighbor in mol.GetAtomWithIdx(carbonyl_carbon).GetNeighbors():
            if neighbor.GetIdx() != oxygen_idx:
                acid_chain_atom = neighbor
                break

        if not acid_chain_atom:
            continue

        # Traverse the alkyl chain from the carbonyl carbon
        chain_length = 0
        current_atom = acid_chain_atom
        visited = set()
        
        while True:
            visited.add(current_atom.GetIdx())
            if current_atom.GetAtomicNum() != 6:
                break
            
            # Check for branching - octanoate is typically straight-chain
            neighbors = [n for n in current_atom.GetNeighbors() 
                        if n.GetIdx() not in visited and n.GetAtomicNum() == 6]
            
            if len(neighbors) != 1:
                break  # branching detected
            
            chain_length += 1
            current_atom = neighbors[0]

        # Octanoate chain has 7 carbons in alkyl chain (total 8 with carbonyl)
        if chain_length >= 7:
            return True, f"Found octanoate chain (length {chain_length+1} carbons)"

    return False, "No octanoate chains found in ester groups"
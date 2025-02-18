"""
Classifies: CHEBI:87658 decanoate ester
"""
"""
Classifies: CHEBI:134741 decanoate ester
"""
from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester is formed by the condensation of decanoic acid (10 carbons including carbonyl)
    with an alcohol/phenol. The acid part must have a straight 10-carbon chain (9 in R group).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # SMARTS pattern for ester group: [OX2][CX3](=[OX1])
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    if not ester_matches:
        return False, "No ester group found"
    
    for match in ester_matches:
        # Match structure: [O]-C(=O)-[O from alcohol]
        carbonyl_idx = match[1]  # Central carbon in ester group
        
        # Find R' carbon (non-oxygen neighbor of carbonyl)
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        r_prime_candidates = []
        for neighbor in carbonyl_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in (match[0], match[2]):
                r_prime_candidates.append(neighbor.GetIdx())
        
        for r_prime_idx in r_prime_candidates:
            # Track chain length and check for straight 9-carbon chain
            current_idx = r_prime_idx
            prev_idx = carbonyl_idx
            chain_length = 1  # Start counting at R' carbon
            valid_chain = True
            
            # Need exactly 9 carbons in straight chain
            for _ in range(8):  # Add 8 more carbons to reach total 9
                current_atom = mol.GetAtomWithIdx(current_idx)
                next_carbons = []
                
                # Find next carbons (non-branching)
                for neighbor in current_atom.GetNeighbors():
                    n_idx = neighbor.GetIdx()
                    if n_idx != prev_idx and neighbor.GetAtomicNum() == 6:
                        next_carbons.append(n_idx)
                
                if len(next_carbons) != 1:
                    valid_chain = False
                    break
                
                prev_idx = current_idx
                current_idx = next_carbons[0]
                chain_length += 1
            
            # Final check: chain length must be exactly 9 and no further carbons
            if valid_chain and chain_length == 9:
                # Verify no additional carbons after 9th
                last_atom = mol.GetAtomWithIdx(current_idx)
                if all(neighbor.GetAtomicNum() != 6 or neighbor.GetIdx() == prev_idx 
                       for neighbor in last_atom.GetNeighbors()):
                    return True, "Straight 10-carbon acid chain (9 in R group)"
    
    return False, "No ester with 10-carbon straight-chain acid"
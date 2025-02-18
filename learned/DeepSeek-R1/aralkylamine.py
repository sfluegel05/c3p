"""
Classifies: CHEBI:18000 aralkylamine
"""
"""
Classifies: aralkylamine
An alkylamine in which the alkyl group is substituted by an aromatic group.
"""
from rdkit import Chem
from rdkit.Chem import Mol

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    An aralkylamine has an amine group connected to an alkyl chain which is substituted by an aromatic group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aralkylamine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Find all amine groups (primary, secondary, tertiary, quaternary; exclude amides)
    # SMARTS: [N] not part of amide, with any number of bonds (X3 or X4 for quaternary)
    amine_pattern = Chem.MolFromSmarts("[N;!$(NC=[O,S,N]);!$(N#*)]")
    amine_matches = mol.GetSubstructMatches(amine_pattern)
    if not amine_matches:
        return False, "No amine group found"
    
    # Check each amine for connection to alkyl chain with aromatic substituent
    for amine_match in amine_matches:
        n_idx = amine_match[0]
        n_atom = mol.GetAtomWithIdx(n_idx)
        
        # Check all carbon neighbors (alkyl chain starts here)
        for neighbor in n_atom.GetNeighbors():
            if neighbor.GetAtomicNum() != 6:
                continue  # Only consider carbon chains
            
            # Perform BFS to traverse all possible alkyl paths from this carbon
            visited = set()
            queue = [(neighbor, n_atom)]  # (current atom, previous atom)
            
            while queue:
                current_atom, prev_atom = queue.pop(0)
                if current_atom.GetIdx() in visited:
                    continue
                visited.add(current_atom.GetIdx())
                
                # Check all neighbors of current atom except previous atom
                for next_atom in current_atom.GetNeighbors():
                    if next_atom.GetIdx() == prev_atom.GetIdx():
                        continue
                    
                    # Check if next atom is part of an aromatic system
                    if next_atom.GetIsAromatic():
                        return True, "Aromatic group substituted on alkyl chain attached to amine"
                    
                    # Check if next atom is in an aromatic ring (even if not marked aromatic)
                    ring_info = mol.GetRingInfo()
                    if next_atom.IsInRing():
                        for ring in ring_info.AtomRings():
                            if next_atom.GetIdx() in ring:
                                # Verify if the ring is aromatic
                                is_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
                                if is_aromatic:
                                    return True, "Aromatic ring substituted on alkyl chain attached to amine"
                
                # Add adjacent carbons to continue traversal
                for next_atom in current_atom.GetNeighbors():
                    if next_atom.GetAtomicNum() == 6 and next_atom.GetIdx() != prev_atom.GetIdx():
                        queue.append((next_atom, current_atom))
    
    return False, "No aromatic substitution on alkyl chain attached to amine"
"""
Classifies: CHEBI:33856 aromatic amino acid
"""
"""
Classifies: aromatic amino acids
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.
    An aromatic amino acid must have an aromatic ring and an amino acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for aromatic ring
    ring_info = mol.GetRingInfo()
    aromatic_rings = []
    for ring in ring_info.AtomRings():
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            aromatic_rings.append(ring)
            
    if not aromatic_rings:
        return False, "No aromatic ring found"

    # Look for alpha-amino acid pattern
    # Matches both primary and N-substituted amino acids
    amino_acid_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0][CX4H1][CX3](=[OX1])[OX2H1]")
    matches = mol.GetSubstructMatches(amino_acid_pattern)
    
    if not matches:
        return False, "No alpha-amino acid group found"

    # For each amino acid group found, verify connection to aromatic ring
    for match in matches:
        n_idx, c_alpha_idx, c_carboxyl_idx, o_idx = match
        
        # Get the alpha carbon atom
        alpha_carbon = mol.GetAtomWithIdx(c_alpha_idx)
        
        # Check all neighbors of alpha carbon and their extended connections
        for neighbor in alpha_carbon.GetNeighbors():
            if neighbor.GetIdx() in (n_idx, c_carboxyl_idx):
                continue
                
            # Check if this branch leads to an aromatic ring
            visited = set()
            stack = [(neighbor, 0)]
            
            while stack:
                current_atom, depth = stack.pop()
                
                if depth > 4:  # Limit search depth
                    continue
                    
                if current_atom.GetIsAromatic():
                    return True, "Contains aromatic ring properly connected to amino acid group"
                    
                visited.add(current_atom.GetIdx())
                
                for next_atom in current_atom.GetNeighbors():
                    if next_atom.GetIdx() not in visited:
                        stack.append((next_atom, depth + 1))
                        
            # Also check if the neighbor itself is part of an aromatic ring
            for ring in aromatic_rings:
                if neighbor.GetIdx() in ring:
                    return True, "Contains aromatic ring properly connected to amino acid group"

    return False, "No proper connection between aromatic ring and amino acid group"
"""
Classifies: CHEBI:29348 fatty amide
"""
"""
Classifies: CHEBI fatty amide (monocarboxylic acid amide derived from a fatty acid)
"""
from rdkit import Chem

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is a monocarboxylic acid amide derived from a fatty acid (aliphatic chain >=4 carbons).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Find all amide groups (CX3=O-NX3)
    amide_pattern = Chem.MolFromSmarts('[CX3](=[OX1])[NX3H0,NX3H1]')
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide group found"

    for match in amide_matches:
        carbonyl_idx = match[0]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Find the R group carbon attached to the carbonyl
        r_carbon = None
        for neighbor in carbonyl_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in {match[1], match[2]}:
                r_carbon = neighbor
                break
        if not r_carbon:
            continue  # No R group (e.g., formamide)

        # Traverse the longest aliphatic chain from R carbon (excluding carbonyl)
        visited = set()
        max_chain = 0
        stack = [(r_carbon, 1)]  # (atom, current_chain_length)
        
        while stack:
            atom, current_length = stack.pop()
            if atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())
            
            if current_length > max_chain:
                max_chain = current_length
            
            # Add adjacent non-aromatic, non-carbonyl carbons
            for neighbor in atom.GetNeighbors():
                if (neighbor.GetAtomicNum() == 6 and 
                    not neighbor.GetIsAromatic() and 
                    neighbor.GetIdx() != carbonyl_idx):
                    stack.append((neighbor, current_length + 1))
        
        # Check if chain length meets minimum for fatty acid (>=4 carbons in R)
        if max_chain >= 4:
            return True, f"Amide with {max_chain}-carbon aliphatic chain"

    return False, "No amide with sufficient aliphatic chain length"
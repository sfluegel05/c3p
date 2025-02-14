"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
"""
Classifies: fatty acid methyl ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    A fatty acid methyl ester is the carboxylic ester obtained by the formal condensation 
    of a fatty acid with methanol.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a fatty acid methyl ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the methyl ester functional group pattern
    # Pattern: Methyl ester group where carbonyl carbon is connected to acyl chain
    methyl_ester_pattern = Chem.MolFromSmarts("COC(=O)[C]")
    ester_matches = mol.GetSubstructMatches(methyl_ester_pattern)
    if not ester_matches:
        return False, "No methyl ester functional group found"
    
    # Analyze each methyl ester group found
    for match in ester_matches:
        # Indices of atoms in the ester group
        methyl_o_idx = match[0]  # O in OCH3
        methyl_c_idx = match[1]  # C in OCH3
        carbonyl_c_idx = match[2]  # Carbonyl C
        acyl_c_idx = match[3]  # First C in acyl chain
        
        # Exclude ester group atoms from acyl chain traversal
        ester_atoms = set(match[:3])

        # Start traversal from the acyl carbon
        acyl_c = mol.GetAtomWithIdx(acyl_c_idx)
        
        # Use BFS to traverse the acyl chain
        visited = set()
        queue = [(acyl_c, None)]  # (atom, parent)
        chain_atoms = set()
        branching_points = 0
        
        # Get ring information
        rings = mol.GetRingInfo().AtomRings()
        cyclic_atoms = set([idx for ring in rings for idx in ring])
        
        while queue:
            atom, parent = queue.pop(0)
            idx = atom.GetIdx()
            if idx in visited or idx in ester_atoms:
                continue
            visited.add(idx)
            chain_atoms.add(idx)
            
            # Check for branching
            neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() != parent and nbr.GetIdx() not in ester_atoms]
            if len(neighbors) > 1:
                branching_points += (len(neighbors) - 1)
            
            # Check for cyclic atoms
            if idx in cyclic_atoms:
                return False, "Acyl chain contains ring structures"
            
            # Check for acceptable elements in acyl chain
            atomic_num = atom.GetAtomicNum()
            # Accept C,H,O,N,S,P, and halogens (F, Cl, Br, I)
            if atomic_num not in (1, 6, 7, 8, 9, 15, 16, 17, 35, 53):
                return False, f"Element {atom.GetSymbol()} not typical in fatty acid chain"
            
            # Add neighbors to queue
            for nbr in neighbors:
                queue.append((nbr, idx))
        
        # Count the number of carbons in the acyl chain
        c_count = sum(1 for idx in chain_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if c_count < 2:
            return False, "Acyl chain too short to be a fatty acid"
        
        # Limit branching points to allow minimal branching
        if branching_points > c_count / 2:
            return False, "Acyl chain too highly branched to be a fatty acid"
        
        # Passed all checks
        return True, "Contains methyl ester group with appropriate fatty acid chain"
    
    return False, "No suitable fatty acid methyl ester pattern found"
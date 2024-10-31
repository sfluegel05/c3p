from rdkit import Chem
from rdkit.Chem import AllChem
from collections import deque

def is_gamma_amino_acid(smiles: str):
    """
    Determines if a molecule is a gamma-amino acid (amino group at gamma position relative to carboxyl).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a gamma-amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Find carboxyl carbons (connected to 2 oxygens, one with double bond)
    carboxyl_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            oxygen_neighbors = []
            double_bond_oxygen = False
            for bond in atom.GetBonds():
                neighbor = bond.GetOtherAtom(atom)
                if neighbor.GetSymbol() == 'O':
                    oxygen_neighbors.append(neighbor)
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        double_bond_oxygen = True
            if len(oxygen_neighbors) == 2 and double_bond_oxygen:
                carboxyl_carbons.append(atom)

    if not carboxyl_carbons:
        return False, "No carboxyl group found"

    # Find amino nitrogens (primary or secondary amines)
    amino_nitrogens = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            h_count = atom.GetTotalNumHs()
            # Only consider primary and secondary amines
            if h_count >= 1:
                amino_nitrogens.append(atom)

    if not amino_nitrogens:
        return False, "No amino group found"

    # For each carboxyl-amino pair, check if they're connected by a valid path
    for carboxyl in carboxyl_carbons:
        carboxyl_idx = carboxyl.GetIdx()
        
        for amino in amino_nitrogens:
            amino_idx = amino.GetIdx()
            
            # Use BFS to find shortest path
            queue = deque([(amino_idx, [amino_idx])])
            visited = {amino_idx}
            valid_path = False
            
            while queue:
                current_idx, path = queue.popleft()
                current = mol.GetAtomWithIdx(current_idx)
                
                if len(path) > 5:  # Path too long
                    continue
                    
                if current_idx == carboxyl_idx and len(path) == 5:
                    # Check if middle atoms are carbons
                    middle_atoms = [mol.GetAtomWithIdx(idx) for idx in path[1:-1]]
                    if all(atom.GetSymbol() == 'C' for atom in middle_atoms):
                        valid_path = True
                        break
                
                for neighbor in current.GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor_idx not in visited:
                        visited.add(neighbor_idx)
                        queue.append((neighbor_idx, path + [neighbor_idx]))

            if valid_path:
                return True, "Valid gamma-amino acid structure found"

    return False, "No valid gamma-amino acid structure found"
# Pr=0.9333333333333333
# Recall=1.0
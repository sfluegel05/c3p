"""
Classifies: CHEBI:35359 carboxamidine
"""
"""
Classifies: CHEBI:51381 carboxamidine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule contains a carboxamidine group (RC(=NR)NR2).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains carboxamidine group, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic carboxamidine pattern: C(=N)N where C is not part of C(=O) or C(=S)
    pattern = Chem.MolFromSmarts("[CX3;!$(C=O);!$(C=S)](=[NX2])[NX3]")
    
    if pattern is None:
        return False, "Error in SMARTS pattern"

    matches = mol.GetSubstructMatches(pattern)
    
    if not matches:
        return False, "No carboxamidine group found"
    
    # For each potential match, do additional validation
    for match in matches:
        carbon_idx = match[0]
        imine_n_idx = match[1]
        amine_n_idx = match[2]
        
        carbon = mol.GetAtomWithIdx(carbon_idx)
        imine_n = mol.GetAtomWithIdx(imine_n_idx)
        amine_n = mol.GetAtomWithIdx(amine_n_idx)
        
        # Exclude guanidines: avoid cases where both nitrogens have additional nitrogen neighbors
        imine_n_neighbors = set(n.GetAtomicNum() for n in imine_n.GetNeighbors())
        amine_n_neighbors = set(n.GetAtomicNum() for n in amine_n.GetNeighbors())
        
        nitrogen_count_imine = sum(1 for x in imine_n_neighbors if x == 7)
        nitrogen_count_amine = sum(1 for x in amine_n_neighbors if x == 7)
        
        if nitrogen_count_imine > 1 or nitrogen_count_amine > 1:
            continue
            
        # Exclude cases where the imine nitrogen is part of N=C-N=C pattern
        # (to avoid certain cyclic imines)
        imine_neighbors = imine_n.GetNeighbors()
        has_additional_imine = False
        for neighbor in imine_neighbors:
            if neighbor.GetIdx() != carbon_idx:  # don't check the carboxamidine carbon
                for next_neighbor in neighbor.GetNeighbors():
                    if next_neighbor.GetAtomicNum() == 7 and next_neighbor.GetIdx() != imine_n_idx:
                        bond = mol.GetBondBetween(neighbor, next_neighbor)
                        if bond.GetBondType() == Chem.BondType.DOUBLE:
                            has_additional_imine = True
                            break
        
        if has_additional_imine:
            continue
            
        # Valid carboxamidine found
        return True, "Contains carboxamidine group"
    
    return False, "No valid carboxamidine group found"
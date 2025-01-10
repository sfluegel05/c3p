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

    # Basic carboxamidine pattern: C(=N)N where:
    # - C is not part of C=O or C=S
    # - The nitrogens can have H or C substituents
    pattern = Chem.MolFromSmarts("[CX3;!$(C=O);!$(C=S);!$(C(=[NH])N=*)](=[NX2])[NX3]")
    
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
        
        # Exclude guanidines by checking if both nitrogens are connected to another nitrogen
        n_neighbors_imine = [n for n in imine_n.GetNeighbors() if n.GetAtomicNum() == 7]
        n_neighbors_amine = [n for n in amine_n.GetNeighbors() if n.GetAtomicNum() == 7]
        
        if len(n_neighbors_imine) > 0 and len(n_neighbors_amine) > 0:
            continue
        
        # Check that carbon has expected number of bonds (should be 3)
        if len(carbon.GetBonds()) != 3:
            continue
            
        # Check that the imine nitrogen is double-bonded to carbon
        bond = mol.GetBondBetweenAtoms(carbon_idx, imine_n_idx)
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
            
        # Check that the amine nitrogen is single-bonded to carbon
        bond = mol.GetBondBetweenAtoms(carbon_idx, amine_n_idx)
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
            
        # Additional pattern to exclude certain cyclic structures where
        # the imine nitrogen is part of another imine group
        cyclic_imine_pattern = Chem.MolFromSmarts("[NX2]=[CX3][NX2]=[CX3]")
        if cyclic_imine_pattern and mol.HasSubstructMatch(cyclic_imine_pattern):
            cyclic_matches = mol.GetSubstructMatches(cyclic_imine_pattern)
            is_part_of_cycle = any(imine_n_idx in match for match in cyclic_matches)
            if is_part_of_cycle:
                continue
        
        # Valid carboxamidine found
        return True, "Contains carboxamidine group (RC(=NR)NR2)"
    
    return False, "No valid carboxamidine group found"
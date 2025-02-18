"""
Classifies: CHEBI:32952 amine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine based on its SMILES string.
    An amine is a compound formally derived from ammonia by replacing one, two, or three hydrogen atoms with hydrocarbyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Preprocess: sanitize and add hydrogens
    mol = Chem.AddHs(mol)
    
    # Check for presence of at least one nitrogen atom
    nitrogens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if not nitrogens:
        return False, "No nitrogen atoms present"
    
    # Patterns to exclude
    nitro_pattern = Chem.MolFromSmarts('[N+](=O)[O-]')  # Nitro group
    amide_pattern = Chem.MolFromSmarts('C(=O)N')        # Amide
    aromatic_n_pattern = Chem.MolFromSmarts('[n]')      # Aromatic nitrogen
    
    for atom in nitrogens:
        # Skip aromatic nitrogens (e.g., pyridine, pyrrole)
        if atom.GetIsAromatic():
            continue
        
        # Check if part of nitro group
        if mol.HasSubstructMatch(nitro_pattern):
            nitro_matches = mol.GetSubstructAtoms(nitro_pattern)
            if atom.GetIdx() in nitro_matches:
                continue
        
        # Check if part of amide group
        if mol.HasSubstructMatch(amide_pattern):
            amide_matches = mol.GetSubstructAtoms(amide_pattern)
            if atom.GetIdx() in amide_matches:
                continue
        
        # Check if adjacent to carbonyl group (like in amides)
        adjacent_to_carbonyl = False
        for bond in atom.GetBonds():
            neighbor = bond.GetOtherAtom(atom)
            if neighbor.GetAtomicNum() == 6:  # Carbon
                for neighbor_bond in neighbor.GetBonds():
                    if neighbor_bond.GetBondType() == Chem.BondType.DOUBLE:
                        other_neighbor = neighbor_bond.GetOtherAtom(neighbor)
                        if other_neighbor.GetAtomicNum() == 8:  # Oxygen
                            adjacent_to_carbonyl = True
                            break
                if adjacent_to_carbonyl:
                    break
        if adjacent_to_carbonyl:
            continue
        
        # Check neighbors: all must be carbon or hydrogen
        carbon_neighbors = []
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 1:  # Hydrogen, allowed
                continue
            elif neighbor.GetAtomicNum() == 6:  # Carbon
                carbon_neighbors.append(neighbor)
        
        # Check number of carbon neighbors (1, 2, or 3)
        num_carbons = len(carbon_neighbors)
        if num_carbons not in {1, 2, 3}:
            continue
        
        # Check total valence (carbons + hydrogens) == 3
        num_h = atom.GetTotalNumHs()
        total_valence = num_carbons + num_h
        if total_valence != 3:
            continue
        
        # All checks passed; it's an amine
        return True, f"Contains {['primary', 'secondary', 'tertiary'][num_carbons-1]} amine group"
    
    # No valid amine groups found
    return False, "No valid amine groups detected"
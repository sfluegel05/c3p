"""
Classifies: CHEBI:26195 polyphenol
"""
from rdkit import Chem

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on its SMILES string.
    A polyphenol is defined as a member of the class of phenols containing 2 or more benzene rings,
    each substituted by at least one hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyphenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find benzene rings
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    benzene_matches = mol.GetSubstructMatches(benzene_pattern)
    
    # Check for hydroxy groups on each benzene ring
    num_benzenes_with_oh = 0
    hydroxy_pattern = Chem.MolFromSmarts("[OH]")
    
    for match in benzene_matches:
        valid_benzene = False
        # Check each atom in the benzene ring
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                # If any adjacent atom is an OH group, mark the benzene
                if mol.HasSubstructMatch(hydroxy_pattern, useChirality=True):
                    valid_benzene = True
                    break
            if valid_benzene:
                break
        
        if valid_benzene:
            num_benzenes_with_oh += 1
            
    if num_benzenes_with_oh >= 2:
        return True, f"Contains {num_benzenes_with_oh} benzene rings each with at least one hydroxy group"
    else:
        return False, f"Found {num_benzenes_with_oh} benzene rings with hydroxy groups, need at least 2"
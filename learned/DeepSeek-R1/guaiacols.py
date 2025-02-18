"""
Classifies: CHEBI:134251 guaiacols
"""
"""
Classifies: CHEBI:28527 guaiacols
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol (phenol with ortho methoxy substituent).
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a guaiacol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Find all phenolic hydroxyl groups (OH directly attached to aromatic carbon)
    phenol_pattern = Chem.MolFromSmarts("[OH]-[c]")
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    if not phenol_matches:
        return False, "No phenolic hydroxyl"
    
    # Check each phenolic OH for adjacent methoxy groups in the same ring
    ortho_methoxy_pattern = Chem.MolFromSmarts("[OH]-[c]1:[c]([OX2]-[CH3]):[c]:[c]:[c]:[c]:1")
    if mol.HasSubstructMatch(ortho_methoxy_pattern):
        return True, "Phenol with ortho methoxy substituent"
    
    # Alternative check for benzene rings with OH and OMe in adjacent positions
    benzene_pattern = Chem.MolFromSmarts("c1:[cH]:[c]([OH]):[c]([OX2]-[CH3]):[cH]:[cH]:1")
    if mol.HasSubstructMatch(benzene_pattern):
        return True, "Phenol with ortho methoxy substituent"
    
    # Check all phenolic positions for adjacent methoxy groups
    for (oh_idx,) in phenol_matches:
        oh_atom = mol.GetAtomWithIdx(oh_idx)
        aromatic_carbon = oh_atom.GetNeighbors()[0]
        
        # Get adjacent atoms in the aromatic system
        for neighbor in aromatic_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 2:
                # Check if oxygen is part of methoxy group
                for o_neighbor in neighbor.GetNeighbors():
                    if o_neighbor.GetAtomicNum() == 6 and o_neighbor.GetTotalNumHs() == 3:
                        return True, "Found methoxy adjacent to phenolic hydroxyl"
    
    return False, "No ortho methoxy relative to phenolic hydroxyl"
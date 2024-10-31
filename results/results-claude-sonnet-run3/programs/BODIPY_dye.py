from rdkit import Chem
from rdkit.Chem import AllChem

def is_BODIPY_dye(smiles: str):
    """
    Determines if a molecule is a BODIPY dye by checking for the presence of the
    4,4-difluoro-4-bora-3a,4a-diaza-s-indacene core structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a BODIPY dye, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # SMARTS pattern for BODIPY core
    # [B-](F)(F) connected to two nitrogens in a specific bicyclic arrangement
    # Note: The pattern explicitly defines the complete core structure
    bodipy_core = Chem.MolFromSmarts('[B-;X4](F)(F)([n,N]1[c,C]2[c,C][c,C][c,C][n,N+]2[c,C][c,C]1)')
    
    # More specific SMARTS pattern that includes the complete core structure
    bodipy_core2 = Chem.MolFromSmarts('[B-;X4](F)(F)([n,N]1[c,C]=,:[c,C][c,C]=,:[c,C][n,N+]2[c,C]=,:[c,C]1)')
    
    if mol.HasSubstructMatch(bodipy_core) or mol.HasSubstructMatch(bodipy_core2):
        # Additional verification of the core structure
        for match in mol.GetSubstructMatches(bodipy_core) + mol.GetSubstructMatches(bodipy_core2):
            boron_idx = match[0]
            boron = mol.GetAtomWithIdx(boron_idx)
            
            # Verify boron has exactly 4 bonds
            if boron.GetDegree() != 4:
                continue
                
            # Verify two fluorines
            f_neighbors = [n for n in boron.GetNeighbors() if n.GetSymbol() == 'F']
            if len(f_neighbors) != 2:
                continue
                
            # Verify two nitrogens
            n_neighbors = [n for n in boron.GetNeighbors() if n.GetSymbol() in ['N', 'n']]
            if len(n_neighbors) != 2:
                continue
            
            # Verify bicyclic system
            ring_info = mol.GetRingInfo()
            if ring_info.NumRings() < 2:
                continue
                
            return True, "Contains BODIPY core structure (4,4-difluoro-4-bora-3a,4a-diaza-s-indacene)"
            
    return False, "Does not contain BODIPY core structure"
# Pr=None
# Recall=None
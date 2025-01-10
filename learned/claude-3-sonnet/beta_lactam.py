"""
Classifies: CHEBI:35627 beta-lactam
"""
"""
Classifies: beta-lactam compounds
A lactam in which the amide bond is contained within a four-membered ring, 
which includes the amide nitrogen and the carbonyl carbon.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a beta-lactam, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Convert to neutral form if possible (to handle charged species)
    mol = Chem.AddHs(mol)
    
    # SMARTS pattern for beta-lactam ring
    # [N;R1]1[C;R1][C;R1](=O)[C;R1]1 represents:
    # - N in a ring (R1) connected to
    # - C in same ring connected to
    # - C(=O) in same ring connected to
    # - C in same ring closing back to N
    beta_lactam_pattern = Chem.MolFromSmarts("[N;R1]1[C;R1][C;R1](=O)[C;R1]1")
    
    # Look for the beta-lactam ring
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"
    
    # Count number of beta-lactam rings
    matches = mol.GetSubstructMatches(beta_lactam_pattern)
    num_rings = len(matches)
    
    # Additional validation - check ring size
    ring_info = mol.GetRingInfo()
    
    # For each match, verify it's part of a 4-membered ring
    valid_rings = 0
    for match in matches:
        ring_atoms = set(match)
        for ring in ring_info.AtomRings():
            if ring_atoms.issubset(set(ring)) and len(ring) == 4:
                valid_rings += 1
                break
    
    if valid_rings == 0:
        return False, "Found amide but not in correct 4-membered ring configuration"
    
    # Additional check for carbonyl
    carbonyl_pattern = Chem.MolFromSmarts("[C;R1](=O)")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "No carbonyl group found in ring"
        
    # Check that nitrogen is part of amide
    amide_n_pattern = Chem.MolFromSmarts("[N;R1]C(=O)")
    if not mol.HasSubstructMatch(amide_n_pattern):
        return False, "Nitrogen not part of amide bond"
    
    return True, f"Contains {valid_rings} beta-lactam ring(s) - 4-membered ring(s) with amide bond"
"""
Classifies: CHEBI:39410 1,2,4-triazines
"""
"""
Classifies: CHEBI:38047 1,2,4-triazine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_2_4_triazines(smiles: str):
    """
    Determines if a molecule contains a 1,2,4-triazine ring structure.
    1,2,4-triazine has nitrogen atoms at positions 1, 2, and 4 of a 6-membered ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains 1,2,4-triazine ring, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for 1,2,4-triazine core
    # [n] represents aromatic nitrogen
    # The numbers ensure the specific 1,2,4 arrangement
    triazine_pattern = Chem.MolFromSmarts("[n]1[n]c[n]cc1")
    
    # Alternative pattern for non-aromatic 1,2,4-triazine
    # This covers cases where the ring is not fully aromatic
    triazine_pattern2 = Chem.MolFromSmarts("N1=NC=NC=C1")
    
    # Additional pattern for tautomeric forms
    triazine_pattern3 = Chem.MolFromSmarts("N1=NC=NN=C1")

    if not (mol.HasSubstructMatch(triazine_pattern) or 
            mol.HasSubstructMatch(triazine_pattern2) or
            mol.HasSubstructMatch(triazine_pattern3)):
        return False, "No 1,2,4-triazine ring found"

    # Additional check to ensure we have a 6-membered ring
    ring_info = mol.GetRingInfo()
    has_six_membered = False
    
    for ring in ring_info.AtomRings():
        if len(ring) == 6:
            ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
            n_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 7)
            if n_count == 3:
                has_six_membered = True
                break
    
    if not has_six_membered:
        return False, "No six-membered ring with exactly 3 nitrogens found"

    # Check aromaticity of the ring
    # Some 1,2,4-triazines might not be aromatic, so this is just informational
    is_aromatic = all(mol.GetAtomWithIdx(i).GetIsAromatic() 
                     for i in mol.GetSubstructMatch(triazine_pattern))
    
    aromatic_status = "aromatic" if is_aromatic else "non-aromatic"
    
    return True, f"Contains {aromatic_status} 1,2,4-triazine ring with nitrogens at positions 1, 2, and 4"
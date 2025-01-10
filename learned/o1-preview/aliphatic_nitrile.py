"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
"""
Classifies: aliphatic nitrile
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile is any nitrile derived from an aliphatic compound.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is an aliphatic nitrile, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude molecules containing metals
    allowed_atomic_nums = {1, 6, 7, 8, 9, 15, 16, 17, 35, 53}  # H, C, N, O, F, P, S, Cl, Br, I
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in allowed_atomic_nums:
            return False, f"Molecule contains disallowed element with atomic number {atomic_num}"
    
    # Check for nitrile group attached to aliphatic carbon
    nitrile_pattern = Chem.MolFromSmarts("[C;X2;H0;+0;!$([C]-[!#6])]#N")
    nitrile_matches = mol.GetSubstructMatches(nitrile_pattern)
    if not nitrile_matches:
        return False, "No nitrile group attached to aliphatic carbon found"
    
    # Ensure the nitrile carbon is not attached to heteroatoms or aromatic atoms
    for match in nitrile_matches:
        nitrile_c_atom = mol.GetAtomWithIdx(match[0])
        if nitrile_c_atom.GetIsAromatic():
            return False, "Nitrile carbon is aromatic, not aliphatic"
        for neighbor in nitrile_c_atom.GetNeighbors():
            if neighbor.GetAtomicNum() != 6:
                return False, "Nitrile carbon attached to heteroatom, not aliphatic carbon"
            if neighbor.GetIsAromatic():
                return False, "Nitrile carbon attached to aromatic carbon, not aliphatic"
    
    # Check for aromatic rings in the molecule
    if mol.GetRingInfo().NumAromaticRings() > 0:
        return False, "Molecule contains aromatic rings, not purely aliphatic"
    
    # Check for aromatic atoms in the molecule
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Molecule contains aromatic atoms, not purely aliphatic"
    
    # If passed all checks, it's an aliphatic nitrile
    return True, "Contains nitrile group attached to aliphatic carbon and molecule is aliphatic"

__metadata__ = {   'chemical_class': {   'name': 'aliphatic nitrile',
                              'definition': 'Any nitrile derived from an aliphatic compound.'},
        'config': {   'llm_model_name': 'lbl/claude-sonnet',
                      'f1_threshold': 0.8}}
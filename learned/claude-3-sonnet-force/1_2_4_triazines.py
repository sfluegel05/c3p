"""
Classifies: CHEBI:39410 1,2,4-triazines
"""
"""
Classifies: CHEBI:25398 1,2,4-triazines
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_2_4_triazine(smiles: str):
    """
    Determines if a molecule is a 1,2,4-triazine based on its SMILES string.
    A 1,2,4-triazine is a compound with a triazine skeleton, where nitrogen atoms
    replace carbon at positions 1, 2, and 4 of the core benzene ring structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1,2,4-triazine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for triazine ring
    triazine_pattern = Chem.MolFromSmarts("c1nncnn1")
    if not mol.HasSubstructMatch(triazine_pattern):
        return False, "No triazine ring found"

    # Check if nitrogen atoms are at positions 1, 2, and 4
    n_atoms = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    n_positions = [i for i, x in enumerate(n_atoms) if x == 7]
    if len(n_positions) != 3:
        return False, f"Wrong number of nitrogen atoms in the ring ({len(n_positions)})"

    ring_info = mol.GetRingInfo()
    triazine_ring = None
    for ring in ring_info.AtomRings():
        if set(ring) == set(n_positions):
            triazine_ring = ring
            break

    if triazine_ring is None:
        return False, "Nitrogen atoms not in ring positions 1, 2, and 4"

    # Additional validation checks
    n_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if n_aromatic_rings == 0:
        return False, "No aromatic rings found"

    n_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable_bonds == 0:
        return False, "No rotatable bonds found (rigid structure)"

    return True, "Contains a 1,2,4-triazine core with nitrogen atoms at positions 1, 2, and 4"
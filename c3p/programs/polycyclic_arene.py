"""
Classifies: CHEBI:33848 polycyclic arene
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene based on its SMILES string.
    A polycyclic arene is a hydrocarbon with at least 3 fused rings and no heteroatoms in the ring system.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polycyclic arene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule is a hydrocarbon
    if not all(atom.GetAtomicNum() in [1, 6] for atom in mol.GetAtoms()):
        return False, "Molecule contains heteroatoms"

    # Check for at least 3 fused rings
    ring_info = mol.GetRingInfo()
    if len(ring_info.AtomRings()) < 3:
        return False, "Molecule does not contain at least 3 fused rings"

    # Check for connected ring system
    ring_system = set()
    for ring in ring_info.AtomRings():
        for atom_idx in ring:
            ring_system.add(atom_idx)
    if len(ring_system) < 10:  # Minimum number of atoms in a polycyclic arene
        return False, "Molecule does not have a connected ring system"

    # Check for aromaticity
    for ring in ring_info.AtomRings():
        ring_mol = Chem.PathToSubmol(mol, ring)
        if not AllChem.IsAromaticRing(ring_mol):
            return False, "Molecule is not aromatic"

    # Check molecular weight and size
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 1000:  # Expected molecular weight range for a polycyclic arene
        return False, "Molecule is outside the expected molecular weight range"

    # Check number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:  # Minimum number of carbon atoms in a polycyclic arene
        return False, "Molecule does not have enough carbon atoms"

    return True, "Molecule is a polycyclic arene"

# Test the function
smiles = "c1ccc-2c(c1)-c1cccc3c4ccccc4cc-2c13"  # benzo[b]fluoranthene
print(is_polycyclic_arene(smiles))
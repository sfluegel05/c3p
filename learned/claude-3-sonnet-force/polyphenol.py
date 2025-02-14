"""
Classifies: CHEBI:26195 polyphenol
"""
"""
Classifies: CHEBI:36328 polyphenol

Polyphenols are defined as members of the class of phenols that contain 2 or more benzene rings,
each of which is substituted by at least one hydroxy group.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aromatic_ring(mol, ring):
    """
    Determines if a ring in a molecule is aromatic.

    Args:
        mol (Mol): RDKit Mol object
        ring (list): List of atom indices in the ring

    Returns:
        bool: True if the ring is aromatic, False otherwise
    """
    # Check if all atoms in the ring are aromatic
    for idx in ring:
        atom = mol.GetAtomWithIdx(idx)
        if not atom.GetIsAromatic():
            return False

    # Check if the ring is planar (all atoms in the same plane)
    ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
    plane = Chem.ComputePlaneThrough(ring_atoms)
    for atom in ring_atoms:
        dist = abs(Chem.ComputeDistanceFrom(atom.GetPoint3D(), plane))
        if dist > 0.1:  # Tolerance of 0.1 Angstrom
            return False

    return True

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on its SMILES string.

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

    # Identify aromatic rings
    aromatic_rings = [ring for ring in Chem.GetSymmSSSR(mol) if is_aromatic_ring(mol, ring)]

    # Check if there are at least 2 benzene rings
    benzene_rings = [ring for ring in aromatic_rings if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring)]
    if len(benzene_rings) < 2:
        return False, "Less than 2 benzene rings"

    # Check for hydroxy groups on benzene rings
    hydroxy_on_benzene = False
    for ring in benzene_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            neighbors = [mol.GetAtomWithIdx(neighbor_idx) for neighbor_idx in atom.GetNeighbors()]
            if any(neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 1 for neighbor in neighbors):
                hydroxy_on_benzene = True
                break
    if not hydroxy_on_benzene:
        return False, "No hydroxy groups on benzene rings"

    return True, "Contains 2 or more benzene rings, with at least one hydroxy group substituted"
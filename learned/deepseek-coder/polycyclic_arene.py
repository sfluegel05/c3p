"""
Classifies: CHEBI:33848 polycyclic arene
"""
"""
Classifies: CHEBI:33853 polycyclic arene
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene based on its SMILES string.
    A polycyclic arene is a polycyclic aromatic hydrocarbon.

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

    # Check if all atoms are aromatic
    for atom in mol.GetAtoms():
        if not atom.GetIsAromatic():
            return False, "Not all atoms are aromatic"

    # Count the number of aromatic rings
    n_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if n_aromatic_rings < 2:
        return False, f"Only {n_aromatic_rings} aromatic ring(s) found, need at least 2"

    # Check if rings are fused (polycyclic)
    ring_info = mol.GetRingInfo()
    if len(ring_info.AtomRings()) < 2:
        return False, "Not enough rings to be polycyclic"

    return True, "Contains multiple fused aromatic rings"
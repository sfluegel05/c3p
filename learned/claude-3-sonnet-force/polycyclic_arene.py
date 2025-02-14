"""
Classifies: CHEBI:33848 polycyclic arene
"""
"""
Classifies: CHEBI:38018 polycyclic arene
A polycyclic aromatic hydrocarbon.
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
    
    # Check aromaticity
    if not mol.GetAromaticForm().IsPossibleAromaticMol():
        return False, "Molecule is not aromatic"
    
    # Count the number of aromatic rings
    aromatic_rings = mol.GetAromaticRings()
    n_aromatic_rings = len(aromatic_rings)
    
    # Polycyclic arenes have at least 2 aromatic rings
    if n_aromatic_rings < 2:
        return False, "Molecule has only one aromatic ring"
    
    # Check if all atoms are carbon or hydrogen
    allowed_atoms = [6, 1]  # Carbon and hydrogen
    atoms = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    if any(atom not in allowed_atoms for atom in atoms):
        return False, "Molecule contains atoms other than carbon and hydrogen"
    
    # Count the number of aromatic bonds
    aromatic_bonds = mol.GetAromaticBondMolWtVector()
    n_aromatic_bonds = len(aromatic_bonds.GetNonzeroElements())
    
    # Polycyclic arenes should have at least 3 aromatic bonds
    if n_aromatic_bonds < 3:
        return False, "Molecule has fewer than 3 aromatic bonds"
    
    # Calculate the ring count
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    
    # Polycyclic arenes should have at least 3 rings
    if ring_count < 3:
        return False, "Molecule has fewer than 3 rings"
    
    return True, "Molecule is a polycyclic aromatic hydrocarbon"
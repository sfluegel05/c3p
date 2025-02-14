"""
Classifies: CHEBI:33848 polycyclic arene
"""
"""
Classifies: CHEBI:33587 polycyclic arene
A polycyclic aromatic hydrocarbon.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Check for aromaticity
    if not mol.GetAromaticForm().IsPossibleAromaticMol():
        return False, "Molecule is not aromatic"

    # Check for polycyclic structure
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 3:
        return False, "Molecule is not polycyclic (fewer than 3 rings)"

    # Check for fused rings
    fused_rings = AllChem.GetSSSR(mol)
    if len(fused_rings) < 2:
        return False, "Molecule does not have fused rings"

    # Check for only carbon and hydrogen atoms
    atoms = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    if any(atom_num not in [1, 6] for atom_num in atoms):
        return False, "Molecule contains atoms other than carbon and hydrogen"

    return True, "Polycyclic aromatic hydrocarbon"
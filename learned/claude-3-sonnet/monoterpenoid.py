"""
Classifies: CHEBI:25409 monoterpenoid
"""
"""
Classifies: CHEBI:24064 monoterpenoid

A monoterpenoid is any terpenoid derived from a monoterpene. The term includes
compounds in which the C10 skeleton of the parent monoterpene has been rearranged
or modified by the removal of one or more skeletal atoms (generally methyl groups).
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count number of atoms
    n_atoms = mol.GetNumAtoms()
    if n_atoms < 10:
        return False, "Too few atoms for monoterpenoid"

    # Check for C10 skeleton
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, "Less than 10 carbon atoms"

    # Check for rearranged or modified monoterpene skeleton
    ring_info = mol.GetRingInfo()
    if ring_info.AtomRings():
        largest_ring_size = max(len(ring) for ring in ring_info.AtomRings())
        if largest_ring_size >= 6:
            return True, "Contains rearranged or modified monoterpene skeleton"

    # Check for characteristic functional groups
    alcohol_pattern = Chem.MolFromSmarts("[OH]")
    ketone_pattern = Chem.MolFromSmarts("C(=O)C")
    ether_pattern = Chem.MolFromSmarts("COC")
    if mol.HasSubstructMatch(alcohol_pattern) or mol.HasSubstructMatch(ketone_pattern) or mol.HasSubstructMatch(ether_pattern):
        return True, "Contains characteristic functional group of monoterpenoids"

    # Check for long hydrocarbon chains
    hydrocarbon_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]")
    if mol.HasSubstructMatch(hydrocarbon_pattern):
        return False, "Likely a sesquiterpenoid or diterpenoid"

    # Default case: not enough evidence to classify
    return None, None
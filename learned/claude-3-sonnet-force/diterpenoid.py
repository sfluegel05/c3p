"""
Classifies: CHEBI:23849 diterpenoid
"""
"""
Classifies: CHEBI:35924 diterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    A diterpenoid is any terpenoid derived from a diterpene, which may have undergone
    rearrangements or modifications of the C20 skeleton, such as the removal of methyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular weight range (more relaxed range)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 500:
        return False, f"Molecular weight {mol_wt:.2f} outside typical diterpenoid range"

    # Look for diterpenoid-like core structure
    core_structure = Chem.MolFromSmarts("[C&r5,r6,r7]")
    if core_structure is None:
        return False, "Invalid SMARTS pattern for diterpenoid core"
    if not mol.HasSubstructMatch(core_structure):
        return False, "No diterpenoid-like core structure found"

    # Look for common functional groups (expanded pattern)
    functional_groups = Chem.MolFromSmarts("[OH,O,C(=O)O,COC,C=C,C=O,C#C]")
    if functional_groups is None:
        return False, "Invalid SMARTS pattern for functional groups"
    if not mol.HasSubstructMatch(functional_groups):
        return False, "No typical diterpenoid functional groups found"

    # Check for bicyclic or tricyclic ring systems
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.BondRings()]
    if not any(size >= 5 and size <= 8 for size in ring_sizes):
        return False, "No bicyclic or tricyclic ring systems found"

    # Count carbon, hydrogen, and oxygen atoms (more relaxed ranges)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if c_count < 15 or c_count > 30:
        return False, f"Carbon count {c_count} outside typical diterpenoid range"
    if h_count < 20 or h_count > 50:
        return False, f"Hydrogen count {h_count} outside typical diterpenoid range"
    if o_count < 0 or o_count > 8:
        return False, f"Oxygen count {o_count} outside typical diterpenoid range"

    return True, "Molecule matches diterpenoid structural features"
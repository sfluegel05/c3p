"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
"""
Classifies: CHEBI:63502 3-oxo-5alpha-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5alpha-steroid based on its SMILES string.
    A 3-oxo-5alpha-steroid is a steroid that has a ketone group at position 3 and
    an alpha configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5alpha-steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid backbone
    steroid_patterns = [
        Chem.MolFromSmarts("[C@]13[C@H]([C@H]([C@@H]2[C@@]([C@@]1(CC[C@H](C2)C)(C)C)(C)CCC(=O)O)C"),
        Chem.MolFromSmarts("[C@]13[C@H]([C@H]([C@@H]2[C@@]([C@@]1(CC[C@H](C2)C)(C)C)(C)CC(=O)O)C"),
        Chem.MolFromSmarts("[C@]13[C@H]([C@H]([C@@H]2[C@@]([C@@]1(CC[C@H](C2)C)(C)C)(C)CCC(=O)C)C"),
        # Add more patterns as needed
    ]
    has_steroid_backbone = any(mol.HasSubstructMatch(pattern) for pattern in steroid_patterns)
    if not has_steroid_backbone:
        return False, "No steroid backbone found"

    # Check for 3-oxo group
    oxo_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds()) == 2]
    has_oxo_group = len(oxo_atoms) > 0

    # Check if oxo group is at position 3
    is_oxo_at_position_3 = False
    for oxo_idx in oxo_atoms:
        oxo_atom = mol.GetAtomWithIdx(oxo_idx)
        neighbor_atoms = [mol.GetAtomWithIdx(nbr_idx) for nbr_idx in oxo_atom.GetNeighbors()]
        if any(nbr_atom.GetAtomicNum() == 6 and nbr_atom.GetDegree() == 4 for nbr_atom in neighbor_atoms):
            is_oxo_at_position_3 = True
            break

    if not has_oxo_group or not is_oxo_at_position_3:
        return False, "No ketone group at position 3"

    # Check for 5alpha configuration
    chiral_centers = Chem.FindMolChiralUnassignedAtoms(mol, includeUnassigned=True)
    has_5alpha_config = any(atom_idx == mol.GetAtomWithIdx(atom_idx).GetNeighbors()[0] - 1 and
                             mol.GetAtomWithIdx(atom_idx).GetHybridization() == Chem.HybridizationType.SP3 and
                             mol.GetAtomWithIdx(atom_idx).GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CW
                             for atom_idx in chiral_centers)

    if not has_5alpha_config:
        return False, "No alpha configuration at position 5"

    # Additional checks for molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 600:
        return False, "Molecular weight outside typical range for steroids"

    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 17 or c_count > 30:
        return False, "Carbon count outside typical range for steroids"

    return True, "Contains 3-oxo and 5alpha configurations in a steroid backbone"
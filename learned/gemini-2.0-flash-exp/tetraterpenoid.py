"""
Classifies: CHEBI:26935 tetraterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    Tetraterpenoids are C40 isoprenoids, often with modifications like methyl loss/rearrangement.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetraterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check carbon count. Must have at least 38 carbons because sometimes methyls are lost
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 38:
        return False, f"Too few carbons: {c_count}. Requires at least 38"

    # Check presence of double bonds. Tetraterpenoids are highly unsaturated.
    double_bond_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) < 5 :
        return False, f"Too few double bonds: {len(double_bond_matches)}. Requires at least 5"

    # Check for methyl groups. Tetraterpenoids are highly methylated
    methyl_pattern = Chem.MolFromSmarts("[CH3X4]")
    methyl_matches = mol.GetSubstructMatches(methyl_pattern)
    if len(methyl_matches) < 6:
        return False, f"Too few methyl groups: {len(methyl_matches)}. Requires at least 6"


    # Check if molecule contains a ring. Cyclic tetraterpenoids are common.
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 1:
        return False, f"No rings found. Many tetraterpenoids contain ring structures"

    # Check molecular weight. Tetraterpenoids tend to be heavy > 500
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low: {mol_wt}. Requires at least 500"

    return True, "Meets tetraterpenoid criteria: C40 isoprenoid structure, multiple double bonds, rings and methyl groups"
"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
"""
Classifies: CHEBI:XXXXX nonclassic icosanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.
    A nonclassic icosanoid is a biologically active signalling molecule made by oxygenation of C20 fatty acids,
    excluding classic icosanoids (leukotrienes and prostanoids).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nonclassic icosanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for C20 backbone (16-24 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (16 <= c_count <= 24):
        return False, f"Not a C20-like molecule (found {c_count} carbons)"

    # Check for oxygenation (at least 1 oxygen atom)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, f"Not enough oxygen atoms (found {o_count})"

    # Check for polyunsaturated pattern (at least 1 double bond)
    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bonds < 1:
        return False, f"Not enough double bonds (found {double_bonds})"

    # Check for functional groups (hydroxyl, epoxy, carboxyl, peroxide, etc.)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    epoxy_pattern = Chem.MolFromSmarts("[C;H1,H2][O][C;H1,H2]")
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    peroxide_pattern = Chem.MolFromSmarts("[OX2][OX2]")
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[C]")

    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    has_epoxy = mol.HasSubstructMatch(epoxy_pattern)
    has_carboxyl = mol.HasSubstructMatch(carboxyl_pattern)
    has_peroxide = mol.HasSubstructMatch(peroxide_pattern)
    has_ketone = mol.HasSubstructMatch(ketone_pattern)

    if not (has_hydroxyl or has_epoxy or has_carboxyl or has_peroxide or has_ketone):
        return False, "Missing characteristic functional groups"

    # Exclude classic icosanoids
    # More specific pattern for leukotrienes (conjugated triene with specific arrangement)
    leukotriene_pattern = Chem.MolFromSmarts("[C]=[C]-[C]=[C]-[C]=[C](-[C])-[C]")
    if mol.HasSubstructMatch(leukotriene_pattern):
        return False, "Contains a conjugated triene system (possible leukotriene)"

    # More specific pattern for prostanoids (cyclopentane with two side chains)
    prostanoid_pattern = Chem.MolFromSmarts("[C]1[C][C][C][C]1(-[C])-[C]")
    if mol.HasSubstructMatch(prostanoid_pattern):
        return False, "Contains a cyclopentane ring with side chains (possible prostanoid)"

    # Check molecular weight (wider range)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 600:
        return False, f"Molecular weight {mol_wt:.2f} is outside the expected range for nonclassic icosanoids"

    return True, "C20-like molecule with oxygenation, double bonds, and functional groups, excluding classic icosanoids"
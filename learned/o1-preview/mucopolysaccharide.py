"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
"""
Classifies: mucopolysaccharide
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide based on its SMILES string.
    A mucopolysaccharide is a polysaccharide composed of repeating units of sugar rings,
    often containing uronic acids and glycosamines, and may be esterified with sulfuric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mucopolysaccharide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Detect 5- or 6-membered rings with at least two heteroatoms (O or N)
    sugar_ring_pattern = Chem.MolFromSmarts("[!R;R5,R6;$(C1~[O,N]~[C,N,O]~[C,N,O]~[C,N,O]~[C,N,O]1)]")
    sugar_rings = mol.GetSubstructMatches(sugar_ring_pattern)
    num_sugar_rings = len(sugar_rings)
    if num_sugar_rings < 2:
        return False, f"Found {num_sugar_rings} sugar-like rings, need at least 2"

    # Detect glycosidic bonds (oxygen bridges between rings)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[C;R]-O-[C;R]")
    glycosidic_bonds = mol.GetSubstructMatches(glycosidic_bond_pattern)
    num_glycosidic_bonds = len(glycosidic_bonds)
    if num_glycosidic_bonds < 1:
        return False, "No glycosidic bonds found between sugar rings"

    # Check for uronic acid units (sugar ring with carboxylic acid substituent)
    uronic_acid_pattern = Chem.MolFromSmarts("[C;R][C;R][C;R][C;R][C;R][O;R][C]=O")
    uronic_acid_units = mol.GetSubstructMatches(uronic_acid_pattern)
    num_uronic_acids = len(uronic_acid_units)
    if num_uronic_acids < 1:
        return False, "No uronic acid units found"

    # Check for glycosamine units (sugar ring with amino group)
    glycosamine_pattern = Chem.MolFromSmarts("[C;R][C;R][C;R][C;R][C;R][O;R][C][N]")
    glycosamine_units = mol.GetSubstructMatches(glycosamine_pattern)
    num_glycosamines = len(glycosamine_units)
    if num_glycosamines < 1:
        return False, "No glycosamine units found"

    # Detect sulfate ester groups (-O-S(=O)(=O)-O-)
    sulfate_ester_pattern = Chem.MolFromSmarts("O-S(=O)(=O)-O")
    sulfate_groups = mol.GetSubstructMatches(sulfate_ester_pattern)
    num_sulfate_groups = len(sulfate_groups)
    # Sulfate groups are common but not mandatory

    # Calculate the total number of rings
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for typical mucopolysaccharide"

    # Check if molecule is a large polymer (due to polysaccharide nature)
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 50:
        return False, "Molecule too small to be a mucopolysaccharide"

    return True, (
        f"Contains {num_sugar_rings} sugar-like rings connected by glycosidic bonds, "
        f"with uronic acid and glycosamine units"
    )

__metadata__ = {
    'chemical_class': {
        'name': 'mucopolysaccharide',
        'definition': 'Any of the group of polysaccharides composed of alternating units from uronic acids and glycosamines, and commonly partially esterified with sulfuric acid.',
    }
}
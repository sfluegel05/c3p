"""
Classifies: CHEBI:18154 polysaccharide
"""
"""
Classifies: CHEBI:36973 polysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    A polysaccharide is a large carbohydrate composed of repeating monosaccharide units 
    linked by glycosidic bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycosidic bonds
    glycosidic_bond_pattern = Chem.MolFromSmarts("[OX2][CX4H][CX4H][OX2]")
    glycosidic_bonds = mol.GetSubstructMatches(glycosidic_bond_pattern)
    if not glycosidic_bonds:
        return False, "No glycosidic bonds found"

    # Check for carbohydrate rings
    carb_ring_pattern = Chem.MolFromSmarts("O[CX4H]1[CX4H]([OX2H])[CX4H]([OX2H])[CX4H]([OX2H])[CX4H]1")
    carb_rings = mol.GetSubstructMatches(carb_ring_pattern)
    if not carb_rings:
        return False, "No carbohydrate rings found"

    # Check for multiple monosaccharide units
    monosaccharide_pattern = Chem.MolFromSmarts("O[CX4H]1[CX4H]([OX2H])[CX4H]([OX2H])[CX4H]([OX2H])[CX4H]1")
    monosaccharide_matches = mol.GetSubstructMatches(monosaccharide_pattern)
    if len(monosaccharide_matches) < 3:
        return False, f"Found {len(monosaccharide_matches)} monosaccharide units, need at least 3"

    # Check molecular weight - polysaccharides typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for polysaccharide"

    return True, "Contains multiple monosaccharide units linked by glycosidic bonds"
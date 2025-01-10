"""
Classifies: CHEBI:26333 prostaglandin
"""
"""
Classifies: prostaglandin compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    Prostaglandins are naturally occurring C20 compounds derived from prostanoic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for cyclopentane ring (5-membered carbon ring)
    cyclopentane_pattern = Chem.MolFromSmarts("[#6]1~[#6]~[#6]~[#6]~[#6]1")
    if not mol.HasSubstructMatch(cyclopentane_pattern):
        return False, "No cyclopentane ring found"

    # Count carbons - should be approximately 20 (allow some variation for derivatives)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15 or c_count > 25:  # Allow some flexibility for derivatives
        return False, f"Carbon count ({c_count}) outside typical range for prostaglandins"

    # Check for presence of oxygen-containing groups (hydroxyls, carbonyls)
    oxygen_pattern = Chem.MolFromSmarts("[OX2H1,O=C]")
    oxygen_matches = mol.GetSubstructMatches(oxygen_pattern)
    if len(oxygen_matches) < 2:  # Prostaglandins typically have multiple oxygen-containing groups
        return False, "Insufficient oxygen-containing groups"

    # Check for carboxylic acid or derivative
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1,OX2]")
    ester_pattern = Chem.MolFromSmarts("[#6][CX3](=[OX1])[OX2][#6]")
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")
    
    has_acid_group = mol.HasSubstructMatch(carboxyl_pattern)
    has_ester = mol.HasSubstructMatch(ester_pattern)
    has_amide = mol.HasSubstructMatch(amide_pattern)
    
    if not (has_acid_group or has_ester or has_amide):
        return False, "No carboxylic acid or derivative found"

    # Check for double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "No carbon-carbon double bonds found"

    # Check for branching from cyclopentane ring
    # We want to see carbons attached to the ring
    ring_with_branches = Chem.MolFromSmarts("[#6]1~[#6](~[#6]~[#6]~[#6]1)~[#6]")
    if not mol.HasSubstructMatch(ring_with_branches):
        return False, "Missing required side chains from cyclopentane ring"

    # Additional check for typical prostaglandin features
    # Look for the common pattern of hydroxyl/carbonyl on the ring
    ring_substituents = Chem.MolFromSmarts("[#6]1~[#6]([OH,=O])~[#6]~[#6]~[#6]1")
    if not mol.HasSubstructMatch(ring_substituents):
        return False, "Missing typical oxygen substitution pattern on ring"

    return True, "Matches prostaglandin structural features: cyclopentane ring, side chains, oxygenated groups, and C20 skeleton"
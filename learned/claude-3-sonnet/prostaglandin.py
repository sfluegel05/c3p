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

    # Count carbons - should be approximately 20 (allow variation for derivatives)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15 or c_count > 25:
        return False, f"Carbon count ({c_count}) outside typical range for prostaglandins"

    # Basic prostaglandin core with cyclopentane ring and alpha chain
    core_pattern = Chem.MolFromSmarts("[CH2,CH]1[CH2,CH][CH2,CH][CH2,CH][CH2,CH]1")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "No cyclopentane ring found"

    # Check for carboxylic acid or derivative at chain end
    acid_patterns = [
        Chem.MolFromSmarts("[CX3](=[OX1])[OX2H]"),  # carboxylic acid
        Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CH2,CH3]"),  # ester
        Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")  # amide
    ]
    
    has_acid_group = any(mol.HasSubstructMatch(patt) for patt in acid_patterns if patt is not None)
    if not has_acid_group:
        return False, "No carboxylic acid or derivative found"

    # Check for oxygen-containing substituents on ring
    ring_oxygen_pattern = Chem.MolFromSmarts("[CH2,CH]1([CH2,CH][CH2,CH]([OH1,=O])[CH2,CH][CH2,CH]1)")
    if not mol.HasSubstructMatch(ring_oxygen_pattern):
        return False, "Missing oxygen substituents on cyclopentane ring"

    # Check for alkenyl (double bond) chains
    alkenyl_pattern = Chem.MolFromSmarts("C=CC")
    if not mol.HasSubstructMatch(alkenyl_pattern):
        return False, "Missing characteristic double bonds"

    # Count total oxygen atoms (typically 4-6 in prostaglandins)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3 or o_count > 8:
        return False, f"Unusual number of oxygen atoms ({o_count}) for prostaglandin"

    # Check for characteristic chain length
    chain_pattern = Chem.MolFromSmarts("CCCC")  # At least 4 carbons in chain
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Missing characteristic carbon chains"

    # Check for hydroxyl groups (common in prostaglandins)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH1]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Missing hydroxyl groups"

    # Additional check for typical prostaglandin features
    # Look for branching from the cyclopentane ring
    branched_ring = Chem.MolFromSmarts("[CH2,CH]1([CH2,CH][CH2,CH][CH2,CH][CH2,CH]1)CC")
    if not mol.HasSubstructMatch(branched_ring):
        return False, "Missing characteristic side chains from cyclopentane ring"

    return True, "Matches prostaglandin structural features: cyclopentane ring with oxygen substituents, characteristic chains and functional groups"
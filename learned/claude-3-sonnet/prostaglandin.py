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
    if c_count < 15 or c_count > 30:
        return False, f"Carbon count ({c_count}) outside typical range for prostaglandins"

    # Basic prostaglandin core with cyclopentane ring
    core_pattern = Chem.MolFromSmarts("[CH2,CH]1[CH2,CH][CH2,CH][CH2,CH][CH2,CH]1")
    if core_pattern is None:
        return None, "Error in core SMARTS pattern"
    if not mol.HasSubstructMatch(core_pattern):
        return False, "No cyclopentane ring found"

    # Look for characteristic chains with double bonds
    chain_pattern = Chem.MolFromSmarts("CC=CC")  # Alkenyl chain
    if chain_pattern is None:
        return None, "Error in chain SMARTS pattern"
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Missing characteristic alkenyl chains"

    # Check for carboxylic acid or derivatives
    acid_pattern = Chem.MolFromSmarts("[$([CX3](=[OX1])[OX2H]),$([CX3](=[OX1])[OX2][CH2,CH3]),$([CX3](=[OX1])[NX3])]")
    if acid_pattern is None:
        return None, "Error in acid SMARTS pattern"
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid or derivative found"

    # Count oxygen atoms (typically 3-8 in prostaglandins and derivatives)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3 or o_count > 8:
        return False, f"Unusual number of oxygen atoms ({o_count}) for prostaglandin"

    # Look for oxygen substituents on or near the ring (hydroxyl, ketone)
    oxygen_pattern = Chem.MolFromSmarts("[$([CH2,CH]1[CH2,CH][CH2,CH]([OH1,=O])[CH2,CH][CH2,CH]1),$([CH2,CH]1[CH2,CH]([OH1,=O])[CH2,CH][CH2,CH][CH2,CH]1)]")
    if oxygen_pattern is None:
        return None, "Error in oxygen SMARTS pattern"
    if not mol.HasSubstructMatch(oxygen_pattern):
        return False, "Missing characteristic oxygen substituents"

    # Check molecular weight - should be in reasonable range for prostaglandins
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 500:
        return False, f"Molecular weight ({mol_wt:.1f}) outside typical range"

    # Check for branching pattern characteristic of prostaglandins
    branch_pattern = Chem.MolFromSmarts("[CH2,CH]1([CH2,CH][CH2,CH][CH2,CH][CH2,CH]1)CC")
    if branch_pattern is None:
        return None, "Error in branch SMARTS pattern"
    if not mol.HasSubstructMatch(branch_pattern):
        return False, "Missing characteristic branching pattern"

    # Additional check for hydroxyl groups (very common in prostaglandins)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH1]")
    if hydroxyl_pattern is None:
        return None, "Error in hydroxyl SMARTS pattern"
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_count < 1:
        return False, "Insufficient hydroxyl groups"

    return True, "Matches prostaglandin structural features: cyclopentane ring, oxygen substituents, characteristic chains and functional groups"
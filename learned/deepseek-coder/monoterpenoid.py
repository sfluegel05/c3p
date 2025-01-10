"""
Classifies: CHEBI:25409 monoterpenoid
"""
"""
Classifies: CHEBI:35189 monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    A monoterpenoid is derived from a monoterpene (C10 skeleton) and may have rearrangements or modifications.

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

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Monoterpenoids typically have 10 carbons, but modifications may expand this range slightly
    if c_count < 8 or c_count > 20:
        return False, f"Carbon count ({c_count}) is outside the expected range for a monoterpenoid"

    # Check for isoprene-like patterns (C5 units) or modified forms
    isoprene_patterns = [
        "[CH3]-[CH2]-[CH]=[CH2]",  # Standard isoprene unit
        "[CH3]-[CH2]-[CH]=[CH]",   # Modified isoprene unit
        "[CH3]-[CH]=[CH2]",        # Simplified isoprene unit
        "[CH3]-[CH]=[CH]",         # Further simplified isoprene unit
        "[CH2]-[CH]=[CH2]",        # Another common pattern
        "[CH2]-[CH]=[CH]",         # Another common pattern
        "[CH3]-[CH2]-[CH2]-[CH]=[CH2]",  # Extended isoprene-like pattern
        "[CH3]-[CH2]-[CH2]-[CH]=[CH]",   # Extended isoprene-like pattern
        "[CH3]-[CH2]-[CH2]-[CH2]-[CH]=[CH2]",  # Extended isoprene-like pattern
        "[CH3]-[CH2]-[CH2]-[CH2]-[CH]=[CH]",   # Extended isoprene-like pattern
    ]
    has_isoprene_pattern = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in isoprene_patterns)
    if not has_isoprene_pattern:
        return False, "No isoprene-like pattern or modified form found"

    # Check for functional groups common in monoterpenoids (e.g., alcohols, ketones, esters, ethers)
    functional_groups = [
        "[OH]",    # Alcohol
        "[C=O]",   # Ketone or aldehyde
        "[O][C=O]",# Ester or carboxylic acid
        "[O][C]",  # Ether
        "[C]=[C]", # Double bond
        "[C#C]",   # Triple bond
        "[C](=O)O" # Carboxylic acid
    ]
    has_functional_group = any(mol.HasSubstructMatch(Chem.MolFromSmarts(fg)) for fg in functional_groups)
    if not has_functional_group:
        return False, "No typical monoterpenoid functional groups found"

    # Check molecular weight (monoterpenoids typically have MW between 130 and 300)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 130 or mol_wt > 300:
        return False, f"Molecular weight ({mol_wt:.2f}) is outside the expected range for a monoterpenoid"

    return True, "Contains a C10 skeleton or modified form with typical monoterpenoid features"
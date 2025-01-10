"""
Classifies: CHEBI:26935 tetraterpenoid
"""
"""
Classifies: CHEBI:25853 tetraterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    Tetraterpenoids are derived from tetraterpenes, which have a C40 skeleton.

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

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Tetraterpenoids typically have around 40 carbon atoms, but allow for modifications
    if c_count < 30 or c_count > 55:
        return False, f"Carbon count ({c_count}) is not consistent with a tetraterpenoid (expected ~40)"

    # Check for extended conjugated systems (characteristic of tetraterpenoids)
    # Look for at least 4 conjugated double bonds in a chain
    conjugated_pattern = Chem.MolFromSmarts("[CX3]=[CX3]~[CX3]=[CX3]~[CX3]=[CX3]")
    conjugated_matches = mol.GetSubstructMatches(conjugated_pattern)
    if len(conjugated_matches) < 1:
        return False, "Insufficient conjugated double bonds for a tetraterpenoid"

    # Check for oxygen-containing functional groups (common in tetraterpenoids)
    oxygen_patterns = [
        "[OX2]",  # Alcohols, ethers
        "[OX1]=[CX3]",  # Carbonyl groups
        "[OX2][CX3](=[OX1])",  # Esters, carboxylic acids
        "[OX2][SX4](=[OX1])(=[OX1])"  # Sulfates
    ]
    oxygen_matches = 0
    for pattern in oxygen_patterns:
        oxygen_matches += len(mol.GetSubstructMatches(Chem.MolFromSmarts(pattern)))
    if oxygen_matches < 1:
        return False, "No oxygen-containing functional groups found"

    # Check molecular weight (tetraterpenoids are typically large molecules)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for a tetraterpenoid"

    # Check for isoprene units (common building blocks of terpenoids)
    isoprene_pattern = Chem.MolFromSmarts("[CX4H2][CX4H]=[CX3H][CX4H2]")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 3:
        return False, "Insufficient isoprene units for a tetraterpenoid"

    return True, "Molecule has a C40-like skeleton with conjugated double bonds, oxygen-containing functional groups, and isoprene units, consistent with a tetraterpenoid"
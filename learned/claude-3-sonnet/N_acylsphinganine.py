"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
"""
Classifies: N-acylsphinganine
A ceramide consisting of sphinganine in which one of the amino hydrogens 
is substituted by a fatty acyl group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an N-acylsphinganine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for amide group (N-acyl)
    amide_pattern = Chem.MolFromSmarts("[NX3H1][CX3](=[OX1])[#6]")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide (N-acyl) group found"
    
    # Check for two hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, f"Found {len(hydroxyl_matches)} hydroxyl groups, need at least 2"
    
    # Check for sphinganine backbone pattern
    # [CH2OH]-[CH]-[NH]-[CH2]n where n is long chain
    sphinganine_pattern = Chem.MolFromSmarts("[CH2X4][OX2H1].[CHX4][OX2H1].[NX3H1]")
    if not mol.HasSubstructMatch(sphinganine_pattern):
        return False, "Missing characteristic sphinganine backbone pattern"
    
    # Count carbons to ensure long enough chains
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18:  # Minimum for sphinganine backbone
        return False, f"Too few carbons ({c_count}) for N-acylsphinganine"
    
    # Check for long alkyl chains using rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Insufficient chain length based on rotatable bonds"
    
    # Verify molecular weight - should be >300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for N-acylsphinganine"
    
    # Count nitrogen atoms - should be exactly one
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count != 1:
        return False, f"Found {n_count} nitrogen atoms, need exactly 1"
    
    # Additional check for primary alcohol (-CH2-OH) and secondary alcohol (-CH-OH)
    primary_alcohol = Chem.MolFromSmarts("[CH2X4][OX2H1]")
    secondary_alcohol = Chem.MolFromSmarts("[CHX4][OX2H1]")
    
    if not (mol.HasSubstructMatch(primary_alcohol) and mol.HasSubstructMatch(secondary_alcohol)):
        return False, "Missing required primary and/or secondary alcohol groups"
        
    # Optional: Check for presence of sugar moieties (common in some N-acylsphinganines)
    sugar_pattern = Chem.MolFromSmarts("[CH1,2]1[OH1,CH2][CH1][CH1][CH1][OH1,CH2]1")
    has_sugar = mol.HasSubstructMatch(sugar_pattern)
    
    base_reason = "Contains sphinganine backbone with N-acyl group and required hydroxyl groups"
    if has_sugar:
        return True, base_reason + " (glycosylated variant)"
    return True, base_reason
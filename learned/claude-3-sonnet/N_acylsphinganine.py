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

    # Key structural features of sphinganine backbone:
    # 1. Long alkyl chain
    # 2. Two hydroxyl groups in specific positions
    # 3. N-acyl group attached to central carbon
    
    # Check for the specific sphinganine backbone pattern
    # [CH2OH]-[CHOH]-[CH]-[NH]-[C=O] connected to long chains
    sphinganine_pattern = Chem.MolFromSmarts(
        "[CH2X4][OX2H1]-[CHX4]([OX2H1])-[CHX4]-[NX3H1]-[CX3](=O)-[#6]"
    )
    if not mol.HasSubstructMatch(sphinganine_pattern):
        return False, "Missing characteristic sphinganine backbone pattern"

    # Verify N-acyl group specifically
    amide_pattern = Chem.MolFromSmarts("[NX3H1][CX3](=[OX1])[CX4]")
    if not mol.GetSubstructMatches(amide_pattern):
        return False, "No proper N-acyl group found"

    # Count carbons to ensure proper chain lengths
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18:
        return False, f"Too few carbons ({c_count}) for N-acylsphinganine"

    # Check for exactly one nitrogen
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count != 1:
        return False, f"Found {n_count} nitrogen atoms, need exactly 1"

    # Count oxygen atoms (minimum 3: two OH groups + one C=O)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, f"Found {o_count} oxygen atoms, need at least 3"

    # Verify presence of both primary and secondary alcohols
    primary_alcohol = Chem.MolFromSmarts("[CH2X4][OX2H1]")
    secondary_alcohol = Chem.MolFromSmarts("[CHX4][OX2H1]")
    
    if not (mol.HasSubstructMatch(primary_alcohol) and mol.HasSubstructMatch(secondary_alcohol)):
        return False, "Missing required primary and/or secondary alcohol groups"

    # Check for minimum chain length using rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
        return False, "Insufficient chain length based on rotatable bonds"

    # Optional: Check for glycosylation
    sugar_pattern = Chem.MolFromSmarts("[CH1,2]1[OH1,CH2][CH1][CH1][CH1][OH1,CH2]1")
    has_sugar = mol.HasSubstructMatch(sugar_pattern)

    # Additional check for fatty acyl chain length
    fatty_chain = Chem.MolFromSmarts("[CX3](=O)[CX4][CX4][CX4][CX4]")
    if not mol.HasSubstructMatch(fatty_chain):
        return False, "Acyl chain too short"

    base_reason = "Contains sphinganine backbone with N-acyl group and required hydroxyl groups"
    if has_sugar:
        return True, base_reason + " (glycosylated variant)"
    return True, base_reason
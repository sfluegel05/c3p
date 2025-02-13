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

    # Core structural requirements:
    # 1. Amide group (N-acyl)
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[#6]")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide (N-acyl) group found"

    # 2. Check for primary alcohol (-CH2-OH) near the amide
    primary_alcohol = Chem.MolFromSmarts("[NX3][CX4][CH2X4][OX2]")
    if not mol.HasSubstructMatch(primary_alcohol):
        return False, "Missing primary alcohol group in correct position"

    # 3. Check for secondary alcohol (-CH(OH)-) near the amide
    secondary_alcohol = Chem.MolFromSmarts("[NX3][CX4][CX4][OX2]")
    if not mol.HasSubstructMatch(secondary_alcohol):
        return False, "Missing secondary alcohol group in correct position"

    # 4. Verify basic sphinganine core structure
    # More flexible pattern allowing for different representations
    sphinganine_core = Chem.MolFromSmarts("[#6][NX3][CX4][CH2X4][OX2]")
    if not mol.HasSubstructMatch(sphinganine_core):
        return False, "Missing characteristic sphinganine core structure"

    # 5. Check for long alkyl chains
    alkyl_chain = Chem.MolFromSmarts("[CX4][CX4][CX4][CX4][CX4]")
    chain_matches = len(mol.GetSubstructMatches(alkyl_chain))
    if chain_matches < 2:
        return False, "Missing required long alkyl chains"

    # Count key atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Basic composition checks
    if c_count < 16:  # Minimum carbon count for shortest N-acylsphinganine
        return False, f"Too few carbons ({c_count}) for N-acylsphinganine"
    if n_count != 1:
        return False, f"Must have exactly 1 nitrogen, found {n_count}"
    if o_count < 3:
        return False, f"Must have at least 3 oxygens, found {o_count}"

    # Check for minimum chain length using rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Insufficient chain length based on rotatable bonds"

    # Check for possible glycosylation
    sugar_pattern = Chem.MolFromSmarts("[CH1,2]1[OH1,CH2][CH1][CH1][CH1][OH1,CH2]1")
    has_sugar = mol.HasSubstructMatch(sugar_pattern)

    base_reason = "Contains sphinganine core with N-acyl group and required hydroxyl groups"
    if has_sugar:
        return True, base_reason + " (glycosylated variant)"
    return True, base_reason
"""
Classifies: CHEBI:50753 isoflavonoid
"""
"""
Classifies: CHEBI:24474 isoflavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is a 1-benzopyran with an aryl substituent at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core pattern for isoflavonoid (3-phenylchromen-4-one core)
    # Specifically matches:
    # - The chromen-4-one core with oxygen at position 1
    # - A phenyl ring at position 3
    # - Allows for various substitutions
    core_pattern = Chem.MolFromSmarts('[#6]1=C(c2[c,C]:[c,C]:[c,C]:[c,C]:[c,C]:2)Oc2c(C1=O):[c,C]:[c,C]:[c,C]:[c,C]:2')
    
    if core_pattern is None:
        return None, "Error in SMARTS pattern"

    # Check for basic isoflavonoid core
    if not mol.HasSubstructMatch(core_pattern):
        return False, "No isoflavonoid core structure found"

    # Additional verification patterns
    
    # Check for chromene oxygen and ketone
    chromene_pattern = Chem.MolFromSmarts('O1c2ccccc2C(=O)C=C1')
    if not mol.HasSubstructMatch(chromene_pattern):
        return False, "Missing required chromene core with ketone"

    # Verify position 3 phenyl substituent
    pos3_phenyl = Chem.MolFromSmarts('C1=COc2ccccc2C1(=O)C1=CC=CC=C1')
    if not mol.HasSubstructMatch(pos3_phenyl):
        return False, "Missing phenyl substituent at position 3"

    # Count ring systems
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Insufficient ring count for isoflavonoid structure"

    # Verify aromatic character
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if aromatic_atoms < 6:
        return False, "Insufficient aromatic character"

    # Additional structural requirements
    required_features = [
        (8, 2, "Insufficient oxygen atoms (need at least 2)"),  # At least 2 oxygens
        (6, 15, "Insufficient carbon atoms for isoflavonoid structure")  # At least 15 carbons
    ]
    
    for atomic_num, min_count, error_msg in required_features:
        atom_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == atomic_num)
        if atom_count < min_count:
            return False, error_msg

    return True, "Contains isoflavonoid core structure (3-phenylchromen-4-one) with required substituents"
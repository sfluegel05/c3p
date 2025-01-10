"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule is an ion or salt and reduce it to the main organic component
    fragments = Chem.GetMolFrags(mol, asMols=True)
    main_mol = max(fragments, key=lambda frag: frag.GetNumAtoms())

    # Check for 6-carbon backbone (hexose based, but allow for derivatives)
    num_carbons = sum(1 for atom in main_mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (6 <= num_carbons <= 8):  # Allow small variations due to possible acetylation, etc.
        return False, f"Number of carbon atoms is {num_carbons}, expected around 6"

    # Check for presence of multiple hydroxyls (account for possible O-acetylation or modifications)
    num_hydroxyls = sum(1 for atom in main_mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1)
    if num_hydroxyls < 3:
        return False, f"Insufficient hydroxyl groups, found {num_hydroxyls} (possible derivatization)"

    # Check for presence of aldehyde or ketone form
    # Aldose pattern: C-C-C-C-C-C(=O)
    # Ketose pattern: C-C-C-C-C(=O)-C
    aldehyde_pattern = Chem.MolFromSmarts("[CX4][CX4][CX4][CX4][CX4][C]=O")
    ketone_pattern = Chem.MolFromSmarts("[CX4][CX4][CX4][CX4][C](=O)[CX4]")
    if not main_mol.HasSubstructMatch(aldehyde_pattern) and not main_mol.HasSubstructMatch(ketone_pattern):
        return False, "Neither aldehyde nor ketone functional groups observed"

    # Adaptation for ring structure - support both cyclic and non-cyclic forms
    ring_info = main_mol.GetRingInfo()
    if not ring_info.NumRings() and not main_mol.HasSubstructMatch(aldehyde_pattern):
        # Hexose in non-cyclic form must have a terminal aldehyde group
        return False, "No ring structures or terminal aldehyde detected, not a typical hexose form"

    return True, "Structure matches criteria for hexose, including variants and derivatives"
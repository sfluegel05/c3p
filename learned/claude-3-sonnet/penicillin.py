"""
Classifies: CHEBI:17334 penicillin
"""
"""
Classifies: CHEBI:17334 penicillin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_penicillin(smiles: str):
    """
    Determines if a molecule is a penicillin based on its SMILES string.
    Penicillins have a characteristic penam core (beta-lactam fused to thiazolidine ring)
    with specific substituents: two methyls at position 2, carboxylate at position 3,
    and carboxamido at position 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a penicillin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for penam core structure (beta-lactam fused to thiazolidine ring)
    # More permissive pattern that matches the basic scaffold
    penam_core = Chem.MolFromSmarts('[#16]1[#6]2[#7][#6](=[#8])[#6][#6]2[#7]1')
    if not mol.HasSubstructMatch(penam_core):
        return False, "No penam core structure found"

    # Check for two methyl groups at position 2
    # More flexible pattern that matches any two methyl groups on the same carbon
    two_methyls = Chem.MolFromSmarts('[CH3][C]([CH3])([#16])[#6]')
    if not mol.HasSubstructMatch(two_methyls):
        return False, "Missing two methyl groups at position 2"

    # Check for carboxylate group at position 3
    # More permissive pattern that matches different forms of carboxylate
    carboxylate = Chem.MolFromSmarts('[#6]([#16])([#6][#7])[#6](=[#8])[#8,#8-]')
    if not mol.HasSubstructMatch(carboxylate):
        return False, "Missing carboxylate group at position 3"

    # Check for carboxamido group at position 6
    # More flexible pattern that matches different types of amide substituents
    carboxamido = Chem.MolFromSmarts('[#6][#7][#6](=[#8])[#6]1[#7][#6](=[#8])')
    if not mol.HasSubstructMatch(carboxamido):
        return False, "Missing carboxamido group at position 6"

    # Basic structural checks
    # Count key atoms to verify overall composition
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    if s_count != 1:
        return False, f"Incorrect number of sulfur atoms (found {s_count}, expected 1)"

    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 2:
        return False, f"Too few nitrogen atoms (found {n_count}, expected at least 2)"

    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:  # At least 3 oxygens (carboxylate + amide)
        return False, f"Too few oxygen atoms (found {o_count}, expected at least 3)"

    # Additional check for bicyclic system size
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Missing required bicyclic system"

    return True, "Contains penam core with correct substituents"
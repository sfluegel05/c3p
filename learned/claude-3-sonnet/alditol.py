"""
Classifies: CHEBI:17522 alditol
"""
"""
Classifies: CHEBI:15972 alditol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    An alditol is an acyclic polyol with formula HOCH2[CH(OH)]nCH2OH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alditol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Reject any molecule with rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        return False, "Contains rings (alditols must be acyclic)"

    # Count atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    other_atoms = sum(1 for atom in mol.GetAtoms() 
                     if atom.GetAtomicNum() not in [1, 6, 8])
    
    if other_atoms > 0:
        return False, "Contains atoms other than C, H, O"
    
    if c_count < 3:
        return False, "Too few carbons for an alditol (minimum 3)"

    # Look for terminal CH2OH groups
    terminal_oh_pattern = Chem.MolFromSmarts("[CH2X4][OX2H1]")
    terminal_oh_matches = len(mol.GetSubstructMatches(terminal_oh_pattern))
    
    if terminal_oh_matches < 2:
        return False, "Missing terminal CH2OH groups"

    # Pattern for internal carbons with hydroxyls: -CH(OH)-
    internal_oh_pattern = Chem.MolFromSmarts("[CHX4]([CX4])[OX2H1]")
    internal_oh_matches = len(mol.GetSubstructMatches(internal_oh_pattern))
    
    # For a valid alditol, each internal carbon should have one hydroxyl
    expected_internal_oh = c_count - 2  # subtract terminal carbons
    if internal_oh_matches < expected_internal_oh:
        return False, "Missing required hydroxyl groups on internal carbons"

    # Check for continuous carbon chain with correct hydroxyl pattern
    # HOCH2-[CHOH]n-CH2OH
    alditol_pattern = Chem.MolFromSmarts("[CH2X4][OX2H1].[CHX4]([OX2H1])[CHX4].[CH2X4][OX2H1]")
    if not mol.HasSubstructMatch(alditol_pattern):
        return False, "Does not match required alditol carbon chain pattern"

    # Count total hydroxyl groups
    oh_groups = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2H1]")))
    
    # For a valid alditol:
    # - Each carbon should have exactly one hydroxyl
    # - Number of hydroxyls should equal number of carbons
    if oh_groups != c_count:
        return False, "Number of hydroxyl groups does not match number of carbons"

    # Additional check for carbonyl groups (should not be present)
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=O)")):
        return False, "Contains carbonyl group (not allowed in alditols)"

    return True, "Valid alditol structure with required hydroxyl pattern"
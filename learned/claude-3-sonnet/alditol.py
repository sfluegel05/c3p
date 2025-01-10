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

    # First check for basic requirements
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 3:
        return False, "Too few carbons for an alditol (minimum 3)"
    if o_count < 3:
        return False, "Too few oxygens for an alditol (minimum 3)"

    # Check for absence of carbonyl groups
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    if mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Contains carbonyl group(s)"

    # Look for terminal CH2OH groups
    terminal_oh_pattern = Chem.MolFromSmarts("[CH2X4][OX2H1]")
    terminal_oh_matches = mol.GetSubstructMatches(terminal_oh_pattern)
    
    # For pure alditols, we need exactly two terminal CH2OH groups
    # However, some derivatives may have modifications, so we'll be lenient here
    if len(terminal_oh_matches) < 1:
        return False, "Missing terminal CH2OH group"

    # Check for chain of carbons with OH groups
    chain_with_oh_pattern = Chem.MolFromSmarts("[CH2X4,CHX4][OX2H1]")
    oh_chain_matches = mol.GetSubstructMatches(chain_with_oh_pattern)
    
    if len(oh_chain_matches) < 3:
        return False, "Insufficient hydroxyl groups in carbon chain"

    # For pure alditols, check if the molecule is mostly linear
    # by comparing the number of branches to the number of atoms
    n_atoms = mol.GetNumAtoms()
    n_branches = 0
    for atom in mol.GetAtoms():
        if atom.GetDegree() > 2:
            n_branches += 1
    
    # Calculate ratio of branching points to total atoms
    branch_ratio = n_branches / n_atoms if n_atoms > 0 else 1
    
    # If it's a pure alditol (not a derivative), it should be mostly linear
    if branch_ratio > 0.5:
        return False, "Too many branches for an alditol"

    # Check for presence of atoms other than C, H, O
    # (allowing for some derivatives that might contain other atoms)
    non_cho_atoms = sum(1 for atom in mol.GetAtoms() 
                       if atom.GetAtomicNum() not in [1, 6, 8])
    
    if non_cho_atoms > 0:
        # This might be an alditol derivative
        return True, "Appears to be an alditol derivative"
    
    # If we've made it here, check the ratio of OH groups to carbons
    # For pure alditols, this should be close to 1
    oh_ratio = o_count / c_count
    
    if 0.8 <= oh_ratio <= 1.2:
        return True, "Pure alditol structure detected"
    else:
        return True, "Possible alditol derivative"
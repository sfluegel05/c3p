"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule is classified as a semisynthetic derivative.
    It attempts to identify natural product cores with significant synthetic modifications.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as semisynthetic, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Diverse natural product core patterns (more coverage)
    core_patterns = [
        Chem.MolFromSmarts("C1=CC=CC=C1"),  # Simple benzene ring
        Chem.MolFromSmarts("C1CCC(CC1)C(=O)O"),  # Cyclohexane carboxylic acid
        Chem.MolFromSmarts("O=C1CC[C@H](N)C2C=CC=CC12"),  # Indole-like
        Chem.MolFromSmarts("C1NC2=C(O1)C=CC=C2"),  # Pyrrolopyridine
        Chem.MolFromSmarts("[#6]1[#6][#6][#6][#6][#6][#6][#6]1"),  # Generic polycyclic structure
    ]
    
    # Broad synthetic modification patterns
    synthetic_patterns = [
        Chem.MolFromSmarts("Cl"),  # Chlorine substitution
        Chem.MolFromSmarts("[C,c][C](=O)[O,N,c]"),  # Amide/ester
        Chem.MolFromSmarts("[C,c]N=[N+]=[N-]"),  # Azide introduction
        Chem.MolFromSmarts("[C,c]S(=O)(=O)[O,N]"),  # Sulfonamide/sulfone
        Chem.MolFromSmarts("[C,c]F"),  # Fluorination
        Chem.MolFromSmarts("[C,c][N+]([O-])=O"),  # Nitro group
        Chem.MolFromSmarts("Br"),  # Bromination
    ]
    
    # Check for natural core patterns
    core_matches = any(mol.HasSubstructMatch(pattern) for pattern in core_patterns)
    
    # Check for synthetic modifications
    synthetic_matches = any(mol.HasSubstructMatch(pattern) for pattern in synthetic_patterns)
    
    # Estimate chemical complexity - an indicator of a synthetic modification
    num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    
    # Evaluate complexity and make a decision
    if core_matches and synthetic_matches:
        reason = "Molecule has natural product-like core with synthetic modifications"
        return True, reason
    
    if num_rings > 3 and num_rotatable_bonds > 5:
        reason = "Molecule shows structural complexity typical of semisynthetic derivatives"
        return True, reason
    
    return False, "No characteristic semisynthetic derivative patterns detected"
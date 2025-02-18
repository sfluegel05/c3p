"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule is classified as a semisynthetic derivative.
    Tries to identify natural product cores with notable synthetic modifications.
    
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

    # Expanded natural product core patterns
    core_patterns = [
        Chem.MolFromSmarts("O[C@H]1[C@@H](O[C@H]2[C@@H](OC1)[C@H](O)[C@@H]1[C@@H]3[C@H](O)[C@@H](O)[C@@H]3C2)C"),  # Complex macrolide
        Chem.MolFromSmarts("C1=CC=CC=C1"),  # Simple aromatic core
        Chem.MolFromSmarts("C1C[C@H]([C@@H]2C[C@H]1O)N2"),  # Beta-lactam structures
        Chem.MolFromSmarts("O=C1C=CC(=O)C2=C1O3[C@@H](C=C2)CC=CC3"),  # Tetracyclic structures
    ]
    
    # Expanded synthetic modification patterns
    synthetic_patterns = [
        Chem.MolFromSmarts("Cl"),  # Halogenation
        Chem.MolFromSmarts("OC(C)=O"),  # Acetylation
        Chem.MolFromSmarts("[C,c](F)"),  # Fluorine additions
        Chem.MolFromSmarts("[N,n]C(=O)OC"),  # Carbamate esters
        Chem.MolFromSmarts("C=N"),  # Nitrile groups
    ]
    
    # Check for natural core patterns
    core_matches = any(mol.HasSubstructMatch(pattern) for pattern in core_patterns)
    
    # Check for synthetic modifications
    synthetic_matches = any(mol.HasSubstructMatch(pattern) for pattern in synthetic_patterns)
    
    # Evaluate structural complexity with natural and synthetic patterns considered
    num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    
    if core_matches and synthetic_matches:
        reason = "Molecule has natural product-like core with synthetic modifications"
        return True, reason
    
    if core_matches and (num_rings >= 4 or num_rotatable_bonds >= 6):
        reason = "Molecule shows structural complexity typical of semisynthetic derivatives"
        return True, reason

    return False, "No characteristic semisynthetic derivative patterns detected"
"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
"""
Classifies: beta-lactam antibiotics
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.
    Beta-lactam antibiotics contain a beta-lactam ring and are used as antibiotics.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam antibiotic, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for beta-lactam ring (4-membered ring with N and C=O)
    # More permissive SMARTS pattern that allows for different substitution patterns
    beta_lactam_pattern = Chem.MolFromSmarts("[NR1]1[CR1][CR1][CR1]1(=O)")
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"

    # Count beta-lactam rings
    beta_lactam_matches = len(mol.GetSubstructMatches(beta_lactam_pattern))
    if beta_lactam_matches > 1:
        return False, "Multiple beta-lactam rings found - unusual for antibiotics"

    # Look for common fused ring systems in beta-lactams with more permissive patterns
    
    # Penicillin-like thiazolidine ring fused to beta-lactam
    penam_pattern = Chem.MolFromSmarts("[NR1]1[CR1]2[CR1][SR1][CR1]2[CR1]1=O")
    
    # Cephalosporin-like 6-membered ring fused to beta-lactam (more permissive)
    cephem_pattern = Chem.MolFromSmarts("[NR1]1[CR1]2[CR1][SR1][CR0,CR1][CR0,CR1]2[CR1]1=O")
    
    # Carbapenem-like 5-membered ring fused to beta-lactam (more permissive)
    carbapenem_pattern = Chem.MolFromSmarts("[NR1]1[CR1]2[CR0,CR1][CR0,CR1][CR1]2[CR1]1=O")
    
    # Monobactam-like structure
    monobactam_pattern = Chem.MolFromSmarts("[NR1]1[CR1][CR1][CR1]1=O")
    
    # Oxacephem-like structure
    oxacephem_pattern = Chem.MolFromSmarts("[NR1]1[CR1]2[CR1][OR1][CR0,CR1][CR0,CR1]2[CR1]1=O")

    # Check for at least one of the common ring systems
    ring_patterns = [penam_pattern, cephem_pattern, carbapenem_pattern, 
                    monobactam_pattern, oxacephem_pattern]
    has_common_ring = any(mol.HasSubstructMatch(pattern) for pattern in ring_patterns)

    # Look for carboxylic acid group (required for activity)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OH,O-]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"

    # Basic size check
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 10:
        return False, "Molecule too small to be a beta-lactam antibiotic"

    # Check molecular weight - beta-lactam antibiotics are typically > 300 Da
    mol_wt = Chem.Descriptors.ExactMolWt(mol)
    if mol_wt < 250:
        return False, "Molecular weight too low for beta-lactam antibiotic"

    return True, "Contains beta-lactam ring and characteristic features of beta-lactam antibiotics"
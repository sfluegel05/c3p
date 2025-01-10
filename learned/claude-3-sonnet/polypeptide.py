"""
Classifies: CHEBI:15841 polypeptide
"""
"""
Classifies: polypeptide (peptide containing ten or more amino acid residues)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide based on its SMILES string.
    A polypeptide is defined as a peptide containing ten or more amino acid residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for peptide bonds pattern [NH]-[C](=O)-[C]
    peptide_pattern = Chem.MolFromSmarts("[NX3H][CX3](=[OX1])[CX4]")
    peptide_matches = mol.GetSubstructMatches(peptide_pattern)
    
    if not peptide_matches:
        return False, "No peptide bonds found"
    
    # Count unique nitrogen atoms involved in peptide bonds
    # This helps avoid counting the same residue multiple times in cyclic peptides
    peptide_nitrogens = set()
    for match in peptide_matches:
        peptide_nitrogens.add(match[0])  # First atom in match is nitrogen
    
    num_residues = len(peptide_nitrogens)
    
    # Additional check for terminal amino group (H2N-CH-)
    terminal_amine_pattern = Chem.MolFromSmarts("[NX3H2][CX4]")
    if mol.HasSubstructMatch(terminal_amine_pattern):
        num_residues += 1
        
    # Check for minimum number of residues (10)
    if num_residues < 10:
        return False, f"Only {num_residues} amino acid residues found, minimum 10 required"
        
    # Additional validation to ensure we're dealing with a peptide structure
    # Look for alpha carbons with attached hydrogens (typical in amino acids)
    alpha_carbon_pattern = Chem.MolFromSmarts("[NX3][CX4H1][CX3](=[OX1])")
    alpha_carbon_matches = mol.GetSubstructMatches(alpha_carbon_pattern)
    
    if len(alpha_carbon_matches) < 5:  # Arbitrary threshold to confirm peptide nature
        return False, "Structure lacks typical amino acid alpha carbon pattern"
    
    # Calculate ratio of peptide bonds to molecule size
    # In peptides, this should be relatively high
    num_atoms = mol.GetNumAtoms()
    peptide_ratio = len(peptide_matches) / num_atoms
    
    if peptide_ratio < 0.1:  # Arbitrary threshold
        return False, "Too few peptide bonds relative to molecule size"
    
    return True, f"Contains {num_residues} amino acid residues connected by peptide bonds"
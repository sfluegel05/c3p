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

    # Multiple patterns for peptide bonds including modified variants
    peptide_patterns = [
        "[NX3H][CX3](=[OX1])[CX4]",  # Standard peptide bond
        "[NX3][CX3](=[OX1])[CX4]",    # N-modified peptide bond
        "[NX3H0][CX3](=[OX1])[CX4]",  # N-substituted peptide bond
    ]
    
    all_peptide_nitrogens = set()
    
    # Find all peptide bonds using multiple patterns
    for pattern in peptide_patterns:
        patt = Chem.MolFromSmarts(pattern)
        matches = mol.GetSubstructMatches(patt)
        for match in matches:
            all_peptide_nitrogens.add(match[0])
    
    # Additional patterns for terminal residues and special cases
    terminal_patterns = [
        "[NX3H2][CX4][CX3](=[OX1])",  # N-terminal
        "[NX3H1][CX4][CX3](=[OX1])",  # Modified N-terminal
        "[CX4][CX3](=[OX1])[OH]",     # C-terminal
        "[CX4][CX3](=[OX1])[O-]",     # C-terminal ionized
    ]
    
    # Check for terminal residues
    terminal_residue_found = False
    for pattern in terminal_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(patt):
            terminal_residue_found = True
            break
            
    # Base count from peptide bonds
    num_residues = len(all_peptide_nitrogens)
    
    # Add one for terminal residue if found
    if terminal_residue_found:
        num_residues += 1
        
    # Check for cyclic peptide patterns
    cyclic_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4]~[NX3][CX3](=[OX1])[CX4]~[NX3][CX3](=[OX1])[CX4]")
    if mol.HasSubstructMatch(cyclic_pattern):
        # Cyclic peptides might need adjustment to residue count
        ring_size = len(mol.GetRingInfo().AtomRings())
        if ring_size > 0:
            num_residues = max(num_residues, len(all_peptide_nitrogens))
    
    # Check for minimum number of residues
    if num_residues < 10:
        return False, f"Only {num_residues} amino acid residues found, minimum 10 required"
    
    # Verify peptide nature by looking for amino acid characteristics
    alpha_carbon_pattern = Chem.MolFromSmarts("[NX3][CX4][CX3](=[OX1])")
    alpha_matches = len(mol.GetSubstructMatches(alpha_carbon_pattern))
    
    if alpha_matches < 5:
        return False, "Structure lacks typical amino acid backbone pattern"
        
    # Success case
    return True, f"Contains {num_residues} amino acid residues connected by peptide bonds"
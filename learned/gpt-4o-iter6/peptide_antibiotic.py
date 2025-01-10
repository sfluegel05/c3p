"""
Classifies: CHEBI:25903 peptide antibiotic
"""
from rdkit import Chem

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.
    Peptide antibiotics are complex peptides with antimicrobial properties, 
    often showing structural motifs like cyclic or branched peptides with 
    unique bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a peptide antibiotic, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for peptide bond pattern: N-C(=O), but in more complex contexts
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)-N")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bond_matches) < 5:
        return False, f"Found {len(peptide_bond_matches)} peptide bonds, less than expected for complex peptide antibiotic"
    
    # Check for unique functional groups typical of peptide antibiotics
    unique_group_patterns = [
        Chem.MolFromSmarts("C1SCCN1"),  # Thiazoline ring
        Chem.MolFromSmarts("OCC2OCNC2"), # Branched ethers
        # Add any other typical peptide antibiotic groups here
    ]
    for pattern in unique_group_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains unique group common in peptide antibiotics"
    
    # Check for cyclic or complex branched structure 
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        return True, "Contains cyclic structure, typical of many peptide antibiotics"

    # Optional: Further domain-specific checks for antibiotic indicators
    
    return False, "No indication of structural motifs typical of peptide antibiotics found"
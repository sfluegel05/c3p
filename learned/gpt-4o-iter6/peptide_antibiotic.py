"""
Classifies: CHEBI:25903 peptide antibiotic
"""
from rdkit import Chem

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.
    Peptide antibiotics are complex peptides with antimicrobial properties, often showing
    structural motifs like cyclic peptides or multiple, consecutive peptide bonds.

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
    
    # Look for peptide bond pattern: N-C(=O)
    peptide_bond_pattern = Chem.MolFromSmarts("N-C(=O)")
    if not mol.HasSubstructMatch(peptide_bond_pattern):
        return False, "No peptide bonds found"

    # Count peptide bonds
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bond_matches) < 5:  # Arbitrary minimum for complexity
        return False, f"Found {len(peptide_bond_matches)} peptide bonds, less than expected for complex peptide antibiotic"
    
    # Check for cyclic structure
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 0:
        return False, "No cyclic structure found, common in peptide antibiotics"
    
    # Optional: Check for diversity of functional groups or unusual amino acids
    # This is hard to define exactly, but a high diversity could be indicative

    return True, "Contains multiple peptide bonds and cyclic structure, typical of peptide antibiotics"

# Testing with example SMILES strings
results = [
    is_peptide_antibiotic("C1CC(=O)NC(C(=O)NC(=O)C2)C1"),  # Simplified structure
    is_peptide_antibiotic("[H]C12NC3=CC=C(C=C3C1(O)CC1N2C(=O)C(NC(=O)C(CC(C)C)NC(=O)C..."),  # Example from the list
]

for result in results:
    print(result)
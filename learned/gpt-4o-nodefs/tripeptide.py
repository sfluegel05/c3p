"""
Classifies: CHEBI:47923 tripeptide
"""
from rdkit import Chem

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide consists of three amino acids linked by peptide bonds.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a tripeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a more comprehensive peptide bond pattern (considering usual common alternatives)
    peptide_bond_pattern = Chem.MolFromSmarts("N[C;R0](=[O;R0])[C;R0]")  # R0 ensures that atoms are not part of a ring
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    
    if len(peptide_bond_matches) != 2:
        return False, f"Found {len(peptide_bond_matches)} peptide bonds, does not match expected tripeptide structure"

    # Check for the presence of a cyclic structure
    ring_info = mol.GetRingInfo()
    
    # If cyclic and contains two peptide bonds, assume it might be a cyclic tripeptide
    if ring_info.NumRings() > 0 and len(peptide_bond_matches) >= 1:
        return True, "Contains cyclic structure with peptide bonds, indicative of a cyclic tripeptide"

    # Validate the number of amino acids (residue-like structures) using a common amino acid backbone pattern
    amino_acid_pattern = Chem.MolFromSmarts("[N;R0][C;R0](=[O;R0])[CX4F,CH3]")  # C-Carbon from residue bay backbone
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)

    if len(amino_acid_matches) != 3:
        return False, f"Found {len(amino_acid_matches)} amino acid residues, does not match expected tripeptide structure"

    return True, "Contains exactly 3 amino acid residues and 2 peptide bonds typical of a linear tripeptide"

# Example usage
# result, reason = is_tripeptide("CC[C@H](C)[C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)N[C@@H](CO)C(O)=O")
# print(result, reason)
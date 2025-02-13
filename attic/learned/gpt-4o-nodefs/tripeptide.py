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

    # Look for peptide bond pattern: N-C(=O)
    peptide_bond_pattern = Chem.MolFromSmarts("N-C(=O)")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    
    # Ensure exactly two peptide bonds for linear tripeptides
    if len(peptide_bond_matches) == 2:
        return True, "Contains exactly 2 peptide bonds typical of a linear tripeptide"

    # Determine if the molecule is cyclic
    kekulized_mol = Chem.Mol(mol)
    contains_ring = kekulized_mol.GetRingInfo().NumRings() > 0

    # If it's cyclic, we may need different logic, but we'll assume cyclic tripeptides for simplicity
    if contains_ring and len(peptide_bond_matches) >= 1:
        return True, f"Contains cyclic structure with peptide bonds, indicative of a cyclic tripeptide"

    # Default to false if none of the expected patterns is properly detected
    return False, f"Found {len(peptide_bond_matches)} peptide bonds, does not match expected tripeptide structure"

# Example usage
# result, reason = is_tripeptide("CC[C@H](C)[C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)N[C@@H](CO)C(O)=O")
# print(result, reason)
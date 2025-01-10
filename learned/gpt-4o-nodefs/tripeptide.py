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
    
    # Find peptide bonds in the molecule
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    num_peptide_bonds = len(peptide_bond_matches)

    # Check for typical amino acid side chains e.g., side chains involved in common amino acids
    # This is a more complex analysis requiring pattern matching for common side chains

    # Check for cyclic structures using SSRS count (number of ring structures in the molecule)
    ssr = Chem.GetSSSR(mol)
    contains_ring = ssr > 0

    # Logic combining these insights:
    # - A linear tripeptide would usually contain two peptide bonds
    # - Cyclic peptides might show different patterns but should also be considered
    if num_peptide_bonds == 2 or contains_ring:
        return True, f"Contains {num_peptide_bonds} peptide bonds typical of a tripeptide structure with possible cyclic structure"

    # Default to false if pattern doesn't align with expected tripeptide characteristics
    return False, f"Found {num_peptide_bonds} peptide bonds, cannot confirm tripeptide without further structure analysis"

# Example usage
# result, reason = is_tripeptide("CC[C@H](C)[C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)N[C@@H](CO)C(O)=O")
# print(result, reason)
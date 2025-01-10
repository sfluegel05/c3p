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

    # Look for the peptide bond pattern: N-C(=O)-C
    peptide_bond_pattern = Chem.MolFromSmarts("N-C(=O)-C")
    
    # Find all the peptide bonds in the molecule
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)

    # Verify if there are exactly three peptide bonds
    if len(peptide_bond_matches) == 3:
        return True, "Contains three peptide bond linkages typical of a tripeptide"
    else:
        return False, f"Found {len(peptide_bond_matches)} peptide bonds, expected exactly 3"

# Example usage
# result, reason = is_tripeptide("CC[C@H](C)[C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)N[C@@H](CO)C(O)=O")
# print(result, reason)
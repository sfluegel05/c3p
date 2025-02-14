"""
Classifies: CHEBI:46761 dipeptide
"""
from rdkit import Chem

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide consists of two amino acid residues linked by a peptide bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dipeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the peptide bond pattern (-C(=O)N-)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    if mol.HasSubstructMatch(peptide_bond_pattern):
        # Locate the peptide bond(s)
        peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
        num_peptide_bonds = len(peptide_bond_matches)

        # A dipeptide generally has one peptide bond
        if num_peptide_bonds == 1:
            return True, "Contains one peptide bond, characteristic of a dipeptide"
        else:
            return False, f"Contains {num_peptide_bonds} peptide bonds, expected exactly 1 for a dipeptide"

    else:
        return False, "No peptide bond found"

# Example usage
# smile = "CC[C@H](C)[C@H](NC(=O)[C@@H](N)CCSC)C(O)=O"  # Example SMILES of Met-Ile
# result = is_dipeptide(smile)
# print(result)
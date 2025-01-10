"""
Classifies: CHEBI:25903 peptide antibiotic
"""
"""
Classifies: peptide antibiotic
"""
from rdkit import Chem

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.
    Peptide antibiotics are characterized by multiple peptide bonds (amide linkages between amino acids).
    This function detects peptide bonds and classifies the molecule accordingly.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a peptide antibiotic, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define peptide bond pattern (amide bond between carbonyl carbon and nitrogen)
    peptide_bond = Chem.MolFromSmarts("C(=O)N")
    if peptide_bond is None:
        return False, "Could not define peptide bond pattern"

    # Find all peptide bonds
    matches = mol.GetSubstructMatches(peptide_bond)
    num_peptide_bonds = len(matches)

    if num_peptide_bonds == 0:
        return False, "No peptide bonds detected"

    # Set a threshold for minimum number of peptide bonds to classify as peptide antibiotic
    if num_peptide_bonds >= 5:
        return True, f"Detected {num_peptide_bonds} peptide bonds"
    else:
        return False, f"Detected {num_peptide_bonds} peptide bonds, which is fewer than required for peptide antibiotic"

__metadata__ = {
    'chemical_class': {
        'name': 'peptide antibiotic',
        'definition': 'A chemically diverse class of peptides that exhibit antimicrobial properties.'
    }
}
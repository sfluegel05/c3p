"""
Classifies: CHEBI:25903 peptide antibiotic
"""
from rdkit import Chem

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.
    A peptide antibiotic is a peptide that exhibits antimicrobial properties.

    Since antimicrobial properties cannot be directly inferred from the SMILES string,
    we will attempt to identify if the molecule is a peptide.

    If the molecule contains multiple peptide bonds, we will consider it as a peptide.
    However, we cannot determine antimicrobial properties from structure alone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a peptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define peptide bond pattern (amide bond between C=O and N)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")

    # Find peptide bonds
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    num_peptide_bonds = len(peptide_bonds)

    if num_peptide_bonds >= 2:
        return True, f"Contains {num_peptide_bonds} peptide bonds indicative of a peptide"
    else:
        return False, "Does not contain sufficient peptide bonds to be a peptide antibiotic"
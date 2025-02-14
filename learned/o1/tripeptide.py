"""
Classifies: CHEBI:47923 tripeptide
"""
"""
Classifies: CHEBI:16699 tripeptide
"""

from rdkit import Chem

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide is an oligopeptide consisting of three amino-acid residues connected by peptide linkages.

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

    # Define a more specific peptide bond pattern
    # Peptide bond: O=C-N-C(alpha), where N is connected to a carbon that is not a carbonyl carbon
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)-N-C")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    num_peptide_bonds = len(peptide_bond_matches)

    # Estimate the number of amino acid residues
    # In linear peptides, Number of residues = Number of peptide bonds + 1
    estimated_residues = num_peptide_bonds + 1

    if estimated_residues != 3:
        return False, f"Estimated {estimated_residues} amino acid residues, need exactly 3 for tripeptide"

    return True, "Molecule is a tripeptide with three amino acid residues connected by peptide bonds"
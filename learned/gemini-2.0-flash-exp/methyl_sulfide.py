"""
Classifies: CHEBI:86315 methyl sulfide
"""
from rdkit import Chem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide based on its SMILES string.
    A methyl sulfide is a sulfide with at least one methyl group attached to the sulfur.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl sulfide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the methyl sulfide SMARTS pattern.  Note that [CX4] means that the carbon has 4 bonds (i.e. is saturated)
    methyl_sulfide_pattern = Chem.MolFromSmarts("[SX2]-[CX4H3]")

    # Check for substructure match.
    if not mol.HasSubstructMatch(methyl_sulfide_pattern):
          return False, "Molecule does not contain a methyl sulfide group"

    # Refine to only aliphatic sulfides
    aliphatic_pattern = Chem.MolFromSmarts("[CX4]-[S]-[CX4]")
    if not mol.HasSubstructMatch(aliphatic_pattern):
        return False, "Molecule is not an aliphatic methyl sulfide"

    return True, "Molecule contains at least one methyl sulfide group and is aliphatic"
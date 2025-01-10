"""
Classifies: CHEBI:86315 methyl sulfide
"""
"""
Classifies: CHEBI:48375 methyl sulfide
"""
from rdkit import Chem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide based on its SMILES string.
    A methyl sulfide is any aliphatic sulfide where at least one of the organyl groups attached to the sulfur is a methyl group.

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

    # Define the methyl sulfide pattern: sulfur bonded to at least one methyl group
    methyl_sulfide_pattern = Chem.MolFromSmarts("[SX2][CX4]")  # S bonded to a carbon with 4 bonds (methyl group)

    # Check if the molecule contains the methyl sulfide pattern
    if mol.HasSubstructMatch(methyl_sulfide_pattern):
        return True, "Contains a sulfur atom bonded to at least one methyl group"
    else:
        return False, "No sulfur atom bonded to a methyl group found"
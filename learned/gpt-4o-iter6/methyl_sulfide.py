"""
Classifies: CHEBI:86315 methyl sulfide
"""
from rdkit import Chem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide based on its SMILES string.
    A methyl sulfide is defined as any aliphatic sulfide in which at least
    one of the organyl groups attached to the sulfur is a methyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl sulfide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Create a SMARTS pattern for sulfur bonded to a methyl group.
    # This allows any other atom or group to be bonded to sulfur.
    methyl_sulfide_pattern = Chem.MolFromSmarts("[CH3][S][*]")

    # Check if the molecule has the substructure pattern
    if mol.HasSubstructMatch(methyl_sulfide_pattern):
        return True, "Contains a sulfur atom bonded to a methyl group (methyl sulfide)"

    return False, "No methyl sulfide pattern found"
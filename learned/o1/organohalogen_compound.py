"""
Classifies: CHEBI:17792 organohalogen compound
"""
"""
Classifies: organohalogen compound
Definition: 'A compound containing at least one carbon-halogen bond (where X is a halogen atom).'
"""
from rdkit import Chem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound based on its SMILES string.
    An organohalogen compound contains at least one carbon-halogen bond (C-X), where X is F, Cl, Br, or I.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organohalogen compound, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for carbon-halogen bond
    c_x_bond = Chem.MolFromSmarts("[#6]-[F,Cl,Br,I]")  # Carbon bonded to halogen

    # Check if molecule has carbon-halogen bond
    if mol.HasSubstructMatch(c_x_bond):
        return True, "Contains at least one carbon-halogen bond"

    # No carbon-halogen bonds found
    return False, "No carbon-halogen bonds found"
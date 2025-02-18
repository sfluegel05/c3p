"""
Classifies: CHEBI:37143 organofluorine compound
"""
"""
Classifies: Organofluorine Compound
Definition: An organofluorine compound is a compound containing at least one carbon-fluorine bond.
"""

from rdkit import Chem

def is_organofluorine_compound(smiles: str):
    """
    Determines if a molecule is an organofluorine compound based on its SMILES string.

    An organofluorine compound is defined as any compound containing at least one carbon-fluorine bond.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule contains a carbon-fluorine bond, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string to a molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for a carbon-fluorine bond.
    # [#6] represents any carbon atom; [F] represents fluorine.
    pattern = Chem.MolFromSmarts("[#6]-[F]")
    
    # Check if the molecule has at least one carbon-fluorine bond.
    if mol.HasSubstructMatch(pattern):
        return True, "Contains at least one carbon-fluorine bond"
    else:
        return False, "No carbon-fluorine bond found"
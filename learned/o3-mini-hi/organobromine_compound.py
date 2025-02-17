"""
Classifies: CHEBI:37141 organobromine compound
"""
"""
Classifies: Organobromine compound
Definition: A compound containing at least one carbon-bromine bond.
"""
from rdkit import Chem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound based on its SMILES string.
    The compound must contain at least one carbon-bromine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organobromine compound, False otherwise
        str: Reason for the classification
    """
    # Parse the SMILES string to an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a carbon-bromine bond.
    # Here, [#6] matches any carbon atom and [Br] matches a bromine atom.
    pattern = Chem.MolFromSmarts("[#6]-[Br]")
    
    # Check if the molecule contains a carbon-bromine bond.
    if not mol.HasSubstructMatch(pattern):
        return False, "No carbon-bromine bond found"
    
    return True, "Contains at least one carbon-bromine bond"
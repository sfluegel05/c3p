"""
Classifies: CHEBI:37141 organobromine compound
"""
from rdkit import Chem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound based on its SMILES string.
    An organobromine compound is defined as a compound containing at least one carbon-bromine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an organobromine compound, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Define the SMARTS pattern for a carbon-bromine single bond
    pattern = Chem.MolFromSmarts("[#6]-Br")  # Carbon atom single-bonded to bromine

    # Check if the molecule contains at least one carbon-bromine bond
    if mol.HasSubstructMatch(pattern):
        return True, "Contains at least one carbon-bromine bond"
    else:
        return False, "Does not contain any carbon-bromine bonds"
"""
Classifies: CHEBI:37141 organobromine compound
"""
from rdkit import Chem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound based on its SMILES string.
    An organobromine compound is defined as an organic compound containing at least one bromine atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an organobromine compound, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule contains at least one bromine atom
    contains_bromine = any(atom.GetAtomicNum() == 35 for atom in mol.GetAtoms())
    if contains_bromine:
        return True, "Contains at least one bromine atom"
    else:
        return False, "Does not contain any bromine atoms"
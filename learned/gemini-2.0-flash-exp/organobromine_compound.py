"""
Classifies: CHEBI:37141 organobromine compound
"""
from rdkit import Chem

def is_organobromine_compound(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is an organobromine compound based on its SMILES string.
    An organobromine compound contains at least one carbon-bromine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple[bool, str]: A tuple containing a boolean indicating if it is an organobromine compound
                          and a string describing the reason.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a C-Br bond using SMARTS
    carbon_bromine_pattern = Chem.MolFromSmarts("[C]Br")
    if mol.HasSubstructMatch(carbon_bromine_pattern):
        return True, "Contains at least one carbon-bromine bond"
    else:
        return False, "Does not contain a carbon-bromine bond"
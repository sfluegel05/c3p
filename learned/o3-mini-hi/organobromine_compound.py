"""
Classifies: CHEBI:37141 organobromine compound
"""
"""
Classifies: Organobromine compound
Definition: A compound containing at least one carbon–bromine bond.
Note: The SMARTS pattern now uses "~" to match any bond type between a carbon and bromine.
"""
from rdkit import Chem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound based on its SMILES string.
    An organobromine compound is defined as a compound containing at least one carbon–bromine bond.
    The SMARTS pattern "[#6]~[Br]" matches any bond (single, aromatic, etc.) between a carbon (atomic number 6)
    and a bromine atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organobromine compound, False otherwise
        str: Reason for the classification
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a carbon–bromine bond.
    # "[#6]" matches any carbon (aromatic or aliphatic) and [Br] matches a bromine atom.
    # The "~" operator allows for any bond type (including aromatic bonds).
    pattern = Chem.MolFromSmarts("[#6]~[Br]")
    
    # Check if there is at least one match of a carbon–bromine bond.
    if not mol.HasSubstructMatch(pattern):
        return False, "No carbon–bromine bond found"
    
    return True, "Contains at least one carbon–bromine bond"
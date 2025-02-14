"""
Classifies: CHEBI:24043 flavones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone based on its SMILES string.
    A flavone has a 2-phenylchromen-4-one skeleton with potential substitutions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavone core substructure using SMARTS
    # [c]1[c](-[c])[o][c](=[O])[c][c]2[c][c][c][c]2[c]1
    flavone_core_smarts = "[c]1[c](-[c])[o][c](=[O])[c][c]2[c][c][c][c]2[c]1"
    flavone_core = Chem.MolFromSmarts(flavone_core_smarts)

    if flavone_core is None:
        return None, "Invalid SMARTS pattern"
    
    # Perform substructure matching
    if mol.HasSubstructMatch(flavone_core):
      return True, "Contains the core flavone structure"
    else:
      return False, "Does not contain the core flavone structure"
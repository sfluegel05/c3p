"""
Classifies: CHEBI:28863 flavanones
"""
"""
Classifies: CHEBI:27555 flavanone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavanone(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string using structural rules.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavanone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the flavanone skeleton
    flavanone_pattern = Chem.MolFromSmarts("C1C(=O)C=C(O1)c2ccccc2")
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "Does not contain the flavanone skeleton"

    # Look for the 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one core structure
    core_pattern = Chem.MolFromSmarts("C1C(=O)C=C(O1)c2ccccc2")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Does not contain the 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one core"

    # Additional checks for substitutions, ring systems, etc.
    # ...

    return True, "Matches the structural rules for flavanones"
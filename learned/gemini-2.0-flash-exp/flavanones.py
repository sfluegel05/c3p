"""
Classifies: CHEBI:28863 flavanones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.
    Flavanones have a 2-aryl-3,4-dihydro-2H-1-benzopyran-4-one skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavanone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for the flavanone core structure
    flavanone_core_smarts = "[C]1[C][C][C]2[O][CH]([CH2][C]2[C]=O)[C]1"
    flavanone_with_aryl_smarts = "[C]1[C][C][C]2[O][CH]([CH2][C]2[C]=O)[C]1[C]3[C][C][C][C][C]3"

    flavanone_core_pattern = Chem.MolFromSmarts(flavanone_core_smarts)
    flavanone_with_aryl_pattern = Chem.MolFromSmarts(flavanone_with_aryl_smarts)

    # Check for the basic flavanone core.
    if not mol.HasSubstructMatch(flavanone_core_pattern):
        return False, "Flavanone core structure not found"
    
    # Check for 2-aryl substitution
    if not mol.HasSubstructMatch(flavanone_with_aryl_pattern):
        return True, "Flavanone core found, but no 2-aryl substitution, likely a flavanone derivative."

    return True, "Molecule is a flavanone"
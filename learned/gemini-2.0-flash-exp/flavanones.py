"""
Classifies: CHEBI:28863 flavanones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.
    Flavanones have a 2-phenyl-3,4-dihydro-2H-1-benzopyran-4-one skeleton and its substituted derivatives.

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

    # Define SMARTS for the flavanone core structure and connected phenyl group
    # Allows for substitutions on the core rings and phenyl group
    flavanone_core_smarts = "[c]1[c][c][c]2[O][CH]([CH2][CH]2[C]=O)[c]1"
    flavanone_with_phenyl_smarts = "[c]1[c][c][c]2[O][CH]([CH2][CH]2[C]=O)[c]1[c]3[c][c][c][c][c]3" # Added aromatic phenyl group

    flavanone_core_pattern = Chem.MolFromSmarts(flavanone_core_smarts)
    flavanone_with_phenyl_pattern = Chem.MolFromSmarts(flavanone_with_phenyl_smarts)

    # Check for the basic flavanone core.
    if not mol.HasSubstructMatch(flavanone_core_pattern):
        return False, "Flavanone core structure not found"
    
    #Check that a phenyl is connected on the 2 position
    if not mol.HasSubstructMatch(flavanone_with_phenyl_pattern):
         return False, "Phenyl group not found at correct position"

    return True, "Molecule is a flavanone"
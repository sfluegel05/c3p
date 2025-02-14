"""
Classifies: CHEBI:71543 rotenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is a rotenoid based on its SMILES string.
    Rotenoids have a tetrahydrochromenochromene skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a rotenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core rotenoid skeleton using SMARTS - a more flexible core
    core_pattern = Chem.MolFromSmarts("[c]12[O][c]3[c]4[c]([c]1)[C](=[O])[c]5[c]([c]([O]3)[c]2)[O][c]45")
    if core_pattern is None:
        return False, "Invalid SMARTS pattern"

    if not mol.HasSubstructMatch(core_pattern):
         return False, "Core tetrahydrochromenochromene skeleton not found"

    return True, "Contains the core tetrahydrochromenochromene skeleton, consistent with a rotenoid"
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

    # Define the core rotenoid skeleton using SMARTS.
    # This pattern aims to capture the fused ring system of rotenoids, relaxing constraints on saturation
    # the 'c' means aromatic, and 'C' means any carbon. ~,- means single or double bonds
    # This is a more relaxed version than previous one. It looks for the 3 fused 6 membered rings, two
    # of them have a bridge oxygen and the carbonyl is always present.
    core_pattern = Chem.MolFromSmarts("[C]1~[C]~[C]~[C](=[O])~[C]2~[c]~[c]~[O]~[c]3~[C]~[C]~[O]~[C]1~[C]3~2")
    
    if core_pattern is None:
        return False, "Invalid SMARTS pattern"


    if not mol.HasSubstructMatch(core_pattern):
         return False, "Core tetrahydrochromenochromene skeleton not found"

    return True, "Contains the core tetrahydrochromenochromene skeleton, consistent with a rotenoid"
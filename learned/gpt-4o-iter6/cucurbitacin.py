"""
Classifies: CHEBI:16219 cucurbitacin
"""
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcNumRings

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are tetracyclic triterpenoids developed by some plants for defense.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a cucurbitacin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for at least four rings to meet the tetracyclic criteria 
    num_rings = CalcNumRings(mol)
    if num_rings < 4:
        return False, "Contains fewer than 4 rings, does not meet tetracyclic criteria"
    
    # Check for functional groups: hydroxyl and carbonyl groups
    hydroxyl_smarts = Chem.MolFromSmarts('[CX4,CX3](O)')
    carbonyl_smarts = Chem.MolFromSmarts('[CX3]=[OX1]')
    
    if not mol.HasSubstructMatch(hydroxyl_smarts):
        return False, "Missing hydroxyl groups in relevant contexts"
    
    if not mol.HasSubstructMatch(carbonyl_smarts):
        return False, "Missing carbonyl groups in the expected context"
    
    # Check for a more generalized cucurbitane-like backbone pattern
    # Rather than precise SMARTS, we'll accommodate variability in the backbone
    cucurbitane_smarts = Chem.MolFromSmarts('C1CCC2C3CC4CCCC(C2)C4C3C1')
    if not mol.HasSubstructMatch(cucurbitane_smarts):
        return False, "Structure does not match a generalized cucurbitane backbone"
    
    return True, "Identified as a cucurbitacin based on structural and functional group analysis"
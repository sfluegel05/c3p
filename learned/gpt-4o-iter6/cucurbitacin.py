"""
Classifies: CHEBI:16219 cucurbitacin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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
    
    # Ensure the molecule has at least 4 rings
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 4:
        return False, "Contains fewer than 4 rings, does not meet tetracyclic criteria"
    
    # Identity ring systems related to cucurbitacins
    cucurbitane_backbone_smarts = [
        Chem.MolFromSmarts('C1[C@H]2CC[C@H]3[C@H]4CC[C@]5([C@]1(CC[C@]23C5)C4)C'),
        Chem.MolFromSmarts('C1=CCC2C3CC4CCC(C2)C4C3CC1'), # More generalized pattern
    ]
    
    if not any(mol.HasSubstructMatch(backbone) for backbone in cucurbitane_backbone_smarts):
        return False, "Structure does not match generalized cucurbitane backbone patterns"
    
    # Check for common cucurbitacin functional group patterns
    hydroxyl_smarts = Chem.MolFromSmarts('[CX4][OX2H]')
    carbonyl_smarts = Chem.MolFromSmarts('[CX3]=[OX1]')
    
    if not mol.HasSubstructMatch(hydroxyl_smarts):
        return False, "Missing hydroxyl groups in relevant contexts"
    
    if not mol.HasSubstructMatch(carbonyl_smarts):
        return False, "Missing carbonyl groups in the expected context"
    
    return True, "Identified as a cucurbitacin based on structural and functional group analysis"
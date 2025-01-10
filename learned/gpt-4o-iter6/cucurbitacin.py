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
    
    # Generalized cucurbitane backbone pattern
    cucurbitane_backbone_smarts = Chem.MolFromSmarts('C1C2CC3(C)C4CCC(C4C5CCC3C5C2C1)C')
    
    if not mol.HasSubstructMatch(cucurbitane_backbone_smarts):
        return False, "Structure does not match generalized cucurbitane backbone pattern"
    
    # Check for common cucurbitacin functional group patterns
    hydroxyl_groups = Chem.MolFromSmarts('[CX4][OX2H]')
    carbonyl_groups = Chem.MolFromSmarts('[CX3]=[OX1]')
    
    if not mol.HasSubstructMatch(hydroxyl_groups):
        return False, "Missing hydroxyl groups in relevant contexts"
    
    if mol.HasSubstructMatch(carbonyl_groups):
        return True, "Identified as a cucurbitacin based on structural and functional group analysis"
    else:
        return False, "Missing necessary carbonyl groups"

    # Consider presence of glycosidic moieties which are common
    glycosidic_moiety_smarts = Chem.MolFromSmarts('O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O') # As an example
    if mol.HasSubstructMatch(glycosidic_moiety_smarts):
        return True, "Contains glycosidic moiety, classified as cucurbitacin"


    return True, "Identified as a cucurbitacin based on detected features"
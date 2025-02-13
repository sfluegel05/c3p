"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.
    The molecule should contain a cyclopentane ring fused to a six-membered oxygen heterocycle,
    indicative of monoterpenoid synthesization from isoprene.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an iridoid monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for cyclopentane ring fused to a six-membered oxygen heterocycle
    iridoid_pattern = Chem.MolFromSmarts("C1CCC2C1O[C@@H]([C@@H]2)")
    secoiridoid_pattern = Chem.MolFromSmarts("O1C=CC(C=C1)")

    is_iridoid = mol.HasSubstructMatch(iridoid_pattern)
    is_secoiridoid = mol.HasSubstructMatch(secoiridoid_pattern)

    if is_iridoid:
        return True, "Contains cyclopentane ring fused to a six-membered oxygen heterocycle, characteristic of an iridoid monoterpenoid"
    
    if is_secoiridoid:
        return True, "Contains open chain structure with bond cleavage, indicative of a secoiridoid monoterpenoid"

    # Further domain-specific heuristics based on secondary structure or modifications can be added

    return False, "Does not match characteristic iridoid monoterpenoid structures"
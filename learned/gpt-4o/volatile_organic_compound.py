"""
Classifies: CHEBI:134179 volatile organic compound
"""
from rdkit import Chem

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound (VOC) based on its SMILES string.
    A VOC is any organic compound with a boiling point <= 250°C under standard atmospheric pressure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a VOC, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check the length of the carbon chain
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count > 20:  # Chains with more than 20 carbons are less likely to be VOCs
        return False, "Carbon chain is likely too long"

    # Check for specific functional groups
    alcohol_pattern = Chem.MolFromSmarts("O")
    ether_pattern = Chem.MolFromSmarts("O-C")
    alkene_pattern = Chem.MolFromSmarts("C=C")
    alkyne_pattern = Chem.MolFromSmarts("C#C")
    halogen_pattern = Chem.MolFromSmarts("[Cl,Br,I,F]")
    
    if mol.HasSubstructMatch(alcohol_pattern):
        return True, "Contains alcohol group"
    if mol.HasSubstructMatch(ether_pattern):
        return True, "Contains ether group"
    if mol.HasSubstructMatch(alkene_pattern):
        return True, "Contains alkene group"
    if mol.HasSubstructMatch(alkyne_pattern):
        return True, "Contains alkyne group"
    if mol.HasSubstructMatch(halogen_pattern):
        return True, "Contains halogen group"

    return False, "Structure does not indicate VOC features"
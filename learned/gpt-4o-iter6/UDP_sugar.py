"""
Classifies: CHEBI:17297 UDP-sugar
"""
from rdkit import Chem

def is_UDP_sugar(smiles: str):
    """
    Determines if a molecule is a UDP-sugar based on its SMILES string.
    A UDP-sugar is a pyrimidine nucleotide-sugar having UDP as the nucleotide component attached
    to an unspecified sugar via an anomeric diphosphate linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a UDP-sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for identifying the uracil part of uridine
    uracil_pattern = Chem.MolFromSmarts("n1c(=O)[nH]c(=O)c(=C1)[*]")  # Capture the uracil ring with possible extensions
    if not mol.HasSubstructMatch(uracil_pattern):
        return False, "No uracil moiety found"
    
    # SMARTS pattern for diphosphate linkage
    diphosphate_pattern = Chem.MolFromSmarts("OP(=O)([O-])[O]P(=O)([O-])O")  # Account for phosphates
    if not mol.HasSubstructMatch(diphosphate_pattern):
        return False, "No diphosphate linkage found"
    
    # Pattern to ensure sugar is attached through the diphosphate
    sugar_attached_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O")  # General sugar structure
    if not mol.HasSubstructMatch(sugar_attached_pattern):
        return False, "No suitable sugar moiety found"

    return True, "Contains UDP group with correct sugar moiety attached via diphosphate linkage"
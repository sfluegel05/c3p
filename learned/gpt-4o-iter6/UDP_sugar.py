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

    # SMARTS pattern for Uridine (nucleoside part of UDP)
    uridine_pattern = Chem.MolFromSmarts("n1c(nc(=O)cc1=O)[*]O[*]O[*]C[*]O[*]")  # Allow variability in linkage
    if not mol.HasSubstructMatch(uridine_pattern):
        return False, "No uridine found"
    
    # Comprehensive SMARTS for the diphosphate and sugar backbone
    # This pattern aims to account for potential charge and connectivity variations
    diphosphate_link_sugar_pattern = Chem.MolFromSmarts("OP(=O)(O)OP(=O)(O)OC[C@H]1O[C@@H](C([C@H](O1)*)*)*")
    
    if not mol.HasSubstructMatch(diphosphate_link_sugar_pattern):
        return False, "No suitable sugar moiety attached via diphosphate linkage"
    
    return True, "Contains UDP group with a suitable sugar moiety attached via diphosphate linkage"
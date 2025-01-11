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
    uridine_pattern = Chem.MolFromSmarts("n1ccc(=O)[nH]c1=O")
    if not mol.HasSubstructMatch(uridine_pattern):
        return False, "No uridine found"

    # SMARTS pattern for diphosphate link (account for potential charge states)
    diphosphate_pattern = Chem.MolFromSmarts("OP(OP(O)(=O)[O-])(=O)OC")
    if not mol.HasSubstructMatch(diphosphate_pattern):
        # Attempt alternative patterns in case of mismatches (neutral states or otherwise)
        diphosphate_pattern_alt = Chem.MolFromSmarts("OP(O)(=O)OP(O)(=O)OC")
        if not mol.HasSubstructMatch(diphosphate_pattern_alt):
            return False, "No diphosphate linkage found"

    # Comprehensive pattern for recognizing UDP-bound sugars
    sugar_moiety_pattern = Chem.MolFromSmarts("OC[C@H]([O-])C([O-])")  # Simplified sugar match, extensible
    diphosphate_and_sugar_pattern = Chem.MolFromSmarts("OP(OP(O)(=O)[O-])(=O)OC[C@@H]1O[C@H]([C@@H]([C@H]1O))")

    # Ensure there's a suitable sugar moiety linked correctly
    if not mol.HasSubstructMatch(diphosphate_and_sugar_pattern):
        return False, "No suitable sugar moiety attached via diphosphate linkage"
    
    return True, "Contains UDP group with a suitable sugar moiety attached via diphosphate linkage"
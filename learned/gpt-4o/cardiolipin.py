"""
Classifies: CHEBI:28494 cardiolipin
"""
from rdkit import Chem

def is_cardiolipin(smiles: str):
    """
    Determines if a molecule is a cardiolipin, defined as a class of phosphatidylglycerols.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a cardiolipin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for two phosphate groups (stronger criteria)
    phosphate_pattern = Chem.MolFromSmarts("OP([O-])([O-])=O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) != 2:
        return False, f"Found {len(phosphate_matches)} phosphate groups, need exactly 2"

    # Check for exactly four ester bonds
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@H]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 4:
        return False, f"Insufficient ester-linked chains, found {len(ester_matches)}, need 4"

    # Check for central glycerol linking
    glycerol_pattern = Chem.MolFromSmarts("C(COP([O-])(=O)O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No central glycerol structure found linking phosphatidic acids"

    # Ensure correct overall connectivity pattern
    cardiolipin_pattern = Chem.MolFromSmarts("C(COP([O-])(=O)O)COCOP([O-])(=O)O")
    if not mol.HasSubstructMatch(cardiolipin_pattern):
        return False, "Structure does not follow cardiolipin connectivity"

    return True, "Molecule structure is consistent with cardiolipin"
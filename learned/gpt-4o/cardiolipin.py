"""
Classifies: CHEBI:28494 cardiolipin
"""
from rdkit import Chem

def is_cardiolipin(smiles: str):
    """
    Determines if a molecule is a cardiolipin.
    
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

    # Check for two phosphate groups
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(O)=O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 2:
        return False, f"Found {len(phosphate_matches)} phosphate groups, need at least 2"

    # Check for ester bonds in long chains
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@H]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 4:
        return False, f"Insufficient ester-linked chains, found {len(ester_matches)}"

    # Check for central glycerol linking phosphates
    linkage_pattern = Chem.MolFromSmarts("C(COP(O)(=O)O)CO")
    linkage_matches = mol.GetSubstructMatches(linkage_pattern)
    if len(linkage_matches) == 0:
        return False, "No central glycerol structure found linking phosphatidic acids"

    return True, "Molecule structure is consistent with cardiolipin"
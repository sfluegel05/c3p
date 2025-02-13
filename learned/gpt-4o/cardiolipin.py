"""
Classifies: CHEBI:28494 cardiolipin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cardiolipin(smiles: str):
    """
    Determines if a molecule is a cardiolipin based on its SMILES string.
    A cardiolipin is a phosphatidylglycerol composed of two molecules of phosphatidic acid covalently linked to a molecule of glycerol.
    
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
    
    # SMARTS pattern for detecting the phosphate group in cardiolipins
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 2:
        return False, f"Found {len(phosphate_matches)} phosphate groups, need at least 2"
    
    # SMARTS pattern for checking long ester-linked chains typically seen in cardiolipins
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@H]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 4:
        return False, f"Insufficient ester-linked chains, found {len(ester_matches)}"

    # Verify a central glycerol connecting two phosphatidic acids through oxygens
    glycerol_pattern = Chem.MolFromSmarts("C(COP(=O)(O)O)CO")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if len(glycerol_matches) == 0:
        return False, "No central glycerol structure found linking two phosphatidyl groups"
    
    return True, "Molecule structure is consistent with cardiolipin"
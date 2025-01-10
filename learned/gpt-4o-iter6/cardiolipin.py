"""
Classifies: CHEBI:28494 cardiolipin
"""
from rdkit import Chem

def is_cardiolipin(smiles: str):
    """
    Determines if a molecule is cardiolipin based on its SMILES string.
    A cardiolipin consists of two phosphatidic acids linked to a glycerol phosphate backbone.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cardiolipin, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern to identify the core glycerol phosphate backbone typical of cardiolipins
    # This should match a central glycerol with phosphate groups
    glycerol_phosphate_pattern = Chem.MolFromSmarts("O[C@@H](COP(=O)(O)O)C(COC(=O)[C@H](O)COP(=O)(O)O)O")

    # SMARTS for phosphatidic acid moiety: phosphate linked with glycerol and esters
    phosphatidic_acid_pattern = Chem.MolFromSmarts("O[C@H](COC(=O)[C@H](O)COP(=O)(O)O)C(=O)O")

    # Check for glycerol phosphate backbone
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "Core glycerol phosphate backbone not found"

    # Find phosphatidic acid moieties
    phosphatidic_acid_matches = mol.GetSubstructMatches(phosphatidic_acid_pattern)
    if len(phosphatidic_acid_matches) < 2:
        return False, f"Found {len(phosphatidic_acid_matches)} phosphatidic acid moieties, need at least 2"

    # Ensure connectivity of glycerol and phosphatidic acids
    glycerol_matches = mol.GetSubstructMatches(glycerol_phosphate_pattern)
    if len(glycerol_matches) != 1:
        return False, "Glycerol backbone isn't connected correctly with phosphatidic acids"

    return True, "Molecule matches cardiolipin structure"
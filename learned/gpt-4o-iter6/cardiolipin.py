"""
Classifies: CHEBI:28494 cardiolipin
"""
from rdkit import Chem

def is_cardiolipin(smiles: str):
    """
    Determines if a molecule is a cardiolipin based on its SMILES string.
    A cardiolipin is defined as a phosphatidylglycerol composed of two molecules 
    of phosphatidic acid covalently linked to a molecule of glycerol.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cardiolipin, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS for identifying the glycerol backbone linked to phosphates
    glycerol_with_phosphate_pattern = Chem.MolFromSmarts("O[C@H](COP(=O)(O)O)[C@@H](O)CO")

    # SMARTS for phosphatidic acid moiety: a phosphate-linked glycerol with ester linkages
    phosphatidic_acid_pattern = Chem.MolFromSmarts("C(=O)O[C@H](COP(O)(=O)O)")

    # Match glycerol backbone connected to phosphates
    if not mol.HasSubstructMatch(glycerol_with_phosphate_pattern):
        return False, "Glycerol backbone with phosphate not found"

    # Ensure there are two phosphatidic acid linkages
    phos_matches = mol.GetSubstructMatches(phosphatidic_acid_pattern)
    if len(phos_matches) < 2:
        return False, f"Found {len(phos_matches)} phosphatidic acid moieties, need at least 2"

    return True, "Molecule matches cardiolipin structure"
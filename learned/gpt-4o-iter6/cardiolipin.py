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

    # SMARTS for glycerol backbone: C-C-C with hydroxyl groups
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")

    # SMARTS for phosphatidic acid moiety: phosphate group linked to glycerol
    phosphatidic_acid_pattern = Chem.MolFromSmarts("COP(O)(O)=O")
    
    # Search for the presence of a glycerol backbone
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Search for at least two phosphatidic acid moieties
    phos_matches = mol.GetSubstructMatches(phosphatidic_acid_pattern)
    if len(phos_matches) < 2:
        return False, f"Found {len(phos_matches)} phosphatidic acid moieties, need at least 2"

    # Further checks to evaluate two sets of glycerols linked through phosphate groups (cardiolipin structure)
    # Assuming enough ester linkages are present based on simplification
    if Chem.MolFromSmarts("COP(=O)(OCCP(OC)=O)OC").HasSubstructMatch(mol):
        return False, "Cardiolipin-like structure not found"

    return True, "Molecule matches cardiolipin structure"
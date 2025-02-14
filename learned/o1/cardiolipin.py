"""
Classifies: CHEBI:28494 cardiolipin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_cardiolipin(smiles: str):
    """
    Determines if a molecule is a cardiolipin based on its SMILES string.
    A cardiolipin consists of a central glycerol connected via phosphodiester bonds
    to two phosphatidic acid units, each with two fatty acid chains.
    
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
    
    # Define SMARTS patterns
    # Central glycerol connected via phosphodiester bonds
    central_glycerol_pattern = Chem.MolFromSmarts("C(COP(O)(=O)O)CO")
    
    # Phosphatidic acid unit (glycerol with ester-linked fatty acids and phosphate)
    phosphatidic_acid_pattern = Chem.MolFromSmarts("C(COP(O)(=O)[O])OC(=O)C")
    
    # Ester-linked fatty acid chains
    fatty_acid_pattern = Chem.MolFromSmarts("C(=O)OCC")
    
    # Check for central glycerol backbone connected via phosphodiester bonds
    glycerol_matches = mol.GetSubstructMatches(central_glycerol_pattern)
    if len(glycerol_matches) == 0:
        return False, "No central glycerol connected via phosphodiester bonds found"
    
    # Check for two phosphatidic acid units
    phosphatidic_acid_matches = mol.GetSubstructMatches(phosphatidic_acid_pattern)
    if len(phosphatidic_acid_matches) < 2:
        return False, f"Found {len(phosphatidic_acid_matches)} phosphatidic acid units, need at least 2"
    
    # Count ester-linked fatty acid chains
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 4:
        return False, f"Found {len(fatty_acid_matches)} ester-linked fatty acid chains, need at least 4"
    
    # Check for phosphate groups
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 2:
        return False, f"Found {len(phosphate_matches)} phosphate groups, need at least 2"
    
    # Additional checks can be implemented to verify the overall structure
    
    return True, "Molecule matches cardiolipin structural features"
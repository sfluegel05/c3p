"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol based on its SMILES string.
    
    A phosphatidylglycerol is defined as having a glycerol backbone esterified by fatty acids,
    with a phosphate group attached to one of the primary hydroxy groups, and another glycerol
    moiety linked via phosphate ester linkage.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a phosphatidylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for identifying a phosphatidylglycerol
    # [O-]P(=O)(OCCO[C@@H](O)CO)OC[C@H](COC=O)OC(=O)C
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](O)CO")
    phosphatidylglycerol_pattern = Chem.MolFromSmarts("O=P(O)(OCCO{}OC={})O{}C{}".format(glycerol_pattern, "C(=O)", "C@", "=O"))
    
    if not mol.HasSubstructMatch(phosphatidylglycerol_pattern):
        return False, "Missing phosphatidylglycerol core structure"
    
    # Verify the presence of a phosphate group (O=P(OH)O)
    phosphate_pattern = Chem.MolFromSmarts("[O-]P(=O)(O)")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Missing phosphate group"
    
    # Verify two esterifiable fatty acid chains are present
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"
    
    # Check for glycerol attachment via phosphate
    glycerol_attachment = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_attachment):
        return False, "Missing glycerol attachment via phosphate"
    
    return True, "Contains a phosphatidylglycerol structure with fatty acid esterified glycerol backbone and phosphate group"
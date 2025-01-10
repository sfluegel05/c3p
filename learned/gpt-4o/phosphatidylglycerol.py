"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol based on its SMILES string.
    
    A phosphatidylglycerol contains a glycerol backbone esterified by fatty acids,
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
    
    # SMARTS patterns for key features in phosphatidylglycerol
    # Glycerol part that should be esterified twice
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Missing glycerol backbone"

    # Phosphate group with glycerol linkage
    phosphoglycerol_pattern = Chem.MolFromSmarts("O=P(O)(OC[C@H](O)CO)OC")
    if not mol.HasSubstructMatch(phosphoglycerol_pattern):
        return False, "Missing phosphoglycerol linkage"
    
    # Ester linkage pattern for fatty acid chains
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"
    
    return True, "Contains phosphatidylglycerol structure with appropriate fatty acid linkages"
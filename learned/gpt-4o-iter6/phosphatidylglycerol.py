"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol based on its SMILES string.

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
    
    # Look for the phosphatidylglycerol pattern
    # Using SMARTS to look for crucial substructures:
    # - Phosphate group attached to the hydroxyl of glycerol
    # - Ester linkages indicating the presence of fatty acid chains

    # Glycerol backbone with a phosphatidyl group pattern (simplified)
    phosphatidylglycerol_pattern = Chem.MolFromSmarts("C([C@@H](COP(O)(=O)O)[O])(CO)O")
    if not mol.HasSubstructMatch(phosphatidylglycerol_pattern):
        return False, "No phosphatidylglycerol backbone found"

    # Check for ester linkages with fatty acid chains
    ester_pattern = Chem.MolFromSmarts("O-C(=O).{2,}")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Insufficient ester linkages, found {len(ester_matches)}"

    # Ensure presence of phosphate group
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Missing phosphate group"

    return True, "Contains glycerol backbone with a phosphatidyl group and ester-linked fatty acid chains"
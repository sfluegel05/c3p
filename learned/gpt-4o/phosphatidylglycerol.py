"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol based on its SMILES string.
    A phosphatidylglycerol is characterized by a glycerol backbone, a phosphatidyl group, and two fatty acid chains.

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

    # Stereochemistry-aware glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts("[C@H](COP(=O)(O)O)[C@@H](CO)O")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Stereochemistry-discrepant or incorrect glycerol backbone"

    # Phosphatidyl group with specificity for PG
    phosphatidyl_pattern = Chem.MolFromSmarts("OP([O-])(=O)O[C@@H](CO)CO")
    if not mol.HasSubstructMatch(phosphatidyl_pattern):
        return False, "The phosphatidyl group is structured differently"

    # Check for two long-chain fatty acid esters
    ester_pattern = Chem.MolFromSmarts("C(=O)OC[C@H](O)CO")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Ester linkages are insufficient, identified {len(ester_matches)}"
    
    # Ensure presence of long hydrophobic chains
    chain_pattern = Chem.MolFromSmarts("CCCCCCC")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if len(chain_matches) < 2:
        return False, f"Long fatty acid chains presence is inadequate, got {len(chain_matches)}"

    return True, "Conforms to phosphatidylglycerol structure with proper glycerol backbone and phosphatidyl group"
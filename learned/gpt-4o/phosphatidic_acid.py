"""
Classifies: CHEBI:16337 phosphatidic acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for glycerol backbone pattern (-OC[CH](OC)-)
    glycerol_pattern = Chem.MolFromSmarts("[O][CH2][CH](O[CH2][O])")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for phosphoric acid group (P(O)(O)=O)
    phosphoric_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")
    if not mol.HasSubstructMatch(phosphoric_pattern):
        return False, "No phosphoric acid group found"
    
    # Look for 2 ester linkages (RCOOR')
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester linkages, need at least 2"

    # Count total number of oxygens to ensure consistency with the structure
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6:
        return False, "Too few oxygens for phosphatidic acid"

    return True, "Contains glycerol backbone, two esterified fatty acid chains, and a phosphoric acid group"
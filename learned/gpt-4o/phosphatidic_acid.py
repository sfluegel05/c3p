"""
Classifies: CHEBI:16337 phosphatidic acid
"""
from rdkit import Chem

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    
    A phosphatidic acid is characterized by a glycerol backbone where two of the hydroxy groups
    are esterified with fatty acids and the third is esterified with phosphoric acid.

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
    
    # Correct glycerol backbone pattern with potential esterified positions
    glycerol_pattern = Chem.MolFromSmarts("C(CO[PH0,PH1]=O)O")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for the presence of phosphoric acid group
    phosphoric_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")
    if not mol.HasSubstructMatch(phosphoric_pattern):
        return False, "No phosphoric acid group found"

    # Check for at least two ester linkages
    ester_pattern = Chem.MolFromSmarts("OC(=O)[C,c]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester linkages, need at least 2"

    # Count total number of oxygens for a complete check
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6:
        return False, "Too few oxygens for phosphatidic acid structure"

    return True, "Contains glycerol backbone, two esterified fatty acid chains, and a phosphoric acid group"
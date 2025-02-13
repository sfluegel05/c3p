"""
Classifies: CHEBI:17855 triglyceride
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triglyceride(smiles: str):
    """
    Determines if a molecule is a triglyceride based on its SMILES string.
    A triglyceride is defined by having a glycerol backbone with three esterified fatty acid chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Enhanced glycerol backbone pattern (with flexible oxygen attachments)
    glycerol_pattern = Chem.MolFromSmarts("[O][CH2][CH](O)[CH2][O]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Ensure three ester groups exist
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 3:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 3"

    # Check for long fatty acid chains
    # Simplified long hydrocarbon chain pattern, allowing for multiple bonds and branches
    fatty_acid_pattern = Chem.MolFromSmarts("C~C~C~C~C~C~C~C")  # A minimum of 8 carbons, general chain
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 3:
        return False, f"Missing sufficient long carbon chains, found {len(fatty_acid_matches)}"

    return True, "Contains glycerol backbone with 3 fatty acid chains attached via ester bonds"
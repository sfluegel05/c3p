"""
Classifies: CHEBI:17389 2-monoglyceride
"""
"""
Classifies: CHEBI:174752 2-monoglyceride
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    A 2-monoglyceride has a single fatty acid esterified at the second position of the glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-monoglyceride, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for exactly one ester group (O-C=O)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, expected 1"
    
    # Check glycerol backbone with ester on C2 and hydroxyls on C1 and C3
    # SMARTS pattern: [CH2]([OH])[CH](O-C=O)[CH2]([OH])
    # More precise pattern to ensure ester is on central carbon with two hydroxyls
    glycerol_pattern = Chem.MolFromSmarts("[CH2]([OH])[CH]([OX2][C](=[OX1])[!O])[CH2]([OH])")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone with ester on C2 and hydroxyls on C1/C3 not found"
    
    return True, "Single ester group on C2 of glycerol with hydroxyls on C1 and C3"
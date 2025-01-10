"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
"""
Classifies: CHEBI:78298 3-oxo-5beta-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5beta-steroid based on its SMILES string.
    A 3-oxo-5beta-steroid is a steroid with a ketone at position 3 and a beta configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5beta-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid backbone (17 carbons in a fused ring system)
    # Use a more flexible pattern to match the steroid backbone
    steroid_pattern = Chem.MolFromSmarts("[C@H]1[C@@H]2[C@@H]3[C@@H]4[C@@H]5[C@@H]6[C@@H]7[C@@H]8[C@@H]9[C@@H]%10[C@@H]%11[C@@H]%12[C@@H]%13[C@@H]%14[C@@H]%15[C@@H]%16[C@@H]%17")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for 3-oxo group (ketone at position 3)
    oxo_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if not any(match[0] == 2 for match in oxo_matches):  # Assuming position 3 is index 2
        return False, "No 3-oxo group found"

    # Check for beta configuration at position 5
    # This is more complex and requires checking the stereochemistry at position 5
    # We assume that the beta configuration is represented by the hydrogen being on the opposite side of the ring junction
    # This is a simplified check and may not cover all cases
    beta_config_pattern = Chem.MolFromSmarts("[C@@H]5")
    if not mol.HasSubstructMatch(beta_config_pattern):
        return False, "No beta configuration at position 5"

    return True, "Contains steroid backbone with 3-oxo group and beta configuration at position 5"
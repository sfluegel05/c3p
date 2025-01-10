"""
Classifies: CHEBI:76575 monoradylglycerol
"""
from rdkit import Chem

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    A monoradylglycerol is a glycerol with one acyl, alkyl, or alk-1-enyl substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("[OH][CX4]([OH])[CX4][OH]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Acyl, alkyl, and alkenyl patterns reflecting varying lengths and stereochemistry
    acyl_pattern = Chem.MolFromSmarts("C(=O)OC(CO)CO")
    alkyl_pattern = Chem.MolFromSmarts("[CH2]OC(CO)CO")
    alkenyl_pattern = Chem.MolFromSmarts("[CH2]=C[CH2]OC(CO)CO")

    # Search for substituents
    acyl_count = len(mol.GetSubstructMatches(acyl_pattern))
    alkyl_count = len(mol.GetSubstructMatches(alkyl_pattern))
    alkenyl_count = len(mol.GetSubstructMatches(alkenyl_pattern))

    # Total distinct substituents should be 1 for a monoradylglycerol
    substituent_count = acyl_count + alkyl_count + alkenyl_count

    if substituent_count == 1:
        return True, "Contains glycerol backbone with one acyl, alkyl, or alk-1-enyl substituent"
    else:
        return False, f"Expected 1 substituent group, found {substituent_count}"
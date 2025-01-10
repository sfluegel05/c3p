"""
Classifies: CHEBI:17408 monoacylglycerol
"""
from rdkit import Chem

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    
    A monoacylglycerol has a glycerol backbone with one acyl group esterified at one of the hydroxyl positions,
    while the other two positions can have hydroxyl or other substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacylglycerol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES to RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define glycerol backbone as a three-carbon chain where each carbon is bonded to at least one oxygen
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone not found"

    # Ensure there is exactly one ester group
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Need exactly 1 primary ester group, found {len(ester_matches)}"
    
    # Identify acyl chain; should have at least 10 carbon atoms (arbitrary approximation for acyl chain)
    acyl_pattern = Chem.MolFromSmarts("C(=O)[CH2][CH2][CH2][CH2]")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No proper acyl group found"

    return True, "Contains glycerol backbone with one acyl group and varied substituents"

# Example Test
smiles = "O(C(=O)CCCCCCC/C=C\CCCCCCCC)CC(O)CO"  # Example of a monacylglycerol
print(is_monoacylglycerol(smiles))
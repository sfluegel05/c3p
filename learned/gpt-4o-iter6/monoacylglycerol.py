"""
Classifies: CHEBI:17408 monoacylglycerol
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

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
        return None, "Invalid SMILES string"
    
    # Define glycerol backbone as C-C-C with at least two oxygens
    glycerol_query = Chem.MolFromSmarts("[C][C](O)[C](O)")
    if not mol.HasSubstructMatch(glycerol_query):
        return False, "Glycerol backbone not found"

    # Define ester group pattern
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    if len(ester_matches) != 1:
        return False, f"Need exactly 1 primary ester group, found {len(ester_matches)}"
    
    # Define Acyl chain pattern
    acyl_pattern = Chem.MolFromSmarts("C(=O)[C,c]")  # Long chain carbon attached pattern
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No proper acyl group found"

    return True, "Contains glycerol backbone with one acyl group and varied substituents"

# Example Test
smiles = "O(C(=O)CCCCCCC/C=C\CCCCCCCC)CC(O)CO"  # Example of a monacylglycerol
print(is_monoacylglycerol(smiles))
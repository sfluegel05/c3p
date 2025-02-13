"""
Classifies: CHEBI:76578 diradylglycerol
"""
from rdkit import Chem

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol is a glycerol backbone with two substituent groups 
    (acyl, alkyl, or alk-1-enyl) attached at any two of the three possible positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for generalized glycerol backbone pattern
    # This pattern looks for -CO- linked carbons with attached oxygens
    glycerol_pattern = Chem.MolFromSmarts("O[C@H](CO)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for the two substituent groups. Each being an ester (-C(=O)O-), alkyl, or alk-1-enyl group.
    # Create patterns for esters, alkyl, and alk-1-enyl 
    ester_pattern = Chem.MolFromSmarts("C(=O)[O]")
    alkyl_enyl_pattern = Chem.MolFromSmarts("C(-C)(-C)(-C)") # Generalized as carbon chain representation
    n_substituents = 0
    
    # Check for ester groups
    n_substituents += len(mol.GetSubstructMatches(ester_pattern))
    # Check for alkyl or alk-1-enyl groups (general hydrocarbon chains)
    n_substituents += len(mol.GetSubstructMatches(alkyl_enyl_pattern))

    # Verify there are two substituent groups (could be any of combinations of above patterns)
    if n_substituents < 2:
        return False, f"Found {n_substituents} substituent groups, need exactly 2"

    return True, "Contains glycerol backbone with two substituent groups"

# Example test case
smiles_example = "C([C@@](COC(CCCCCC/C=C\\C/C=C\\C/C=C\\CCCCC)=O)(OC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)[H])O"
print(is_diradylglycerol(smiles_example))
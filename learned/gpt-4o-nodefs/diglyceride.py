"""
Classifies: CHEBI:18035 diglyceride
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_diglyceride(smiles: str):
    """
    Determines if a molecule is a diglyceride based on its SMILES string.
    A diglyceride is a glycerol backbone with two fatty acid chains attached via ester bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with hydroxyl groups)
    glycerol_pattern = Chem.MolFromSmarts("[CX3](OC)[CX3](OC)[CX3](O)")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with two hydroxyl groups found"
        
    # Look for 2 ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[C](=O)[OX2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Check for at least two long carbon chains (fatty acids)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, "Too few long carbon chains for fatty acids"

    return True, "Contains glycerol backbone with two fatty acid chains attached via ester bonds"

## Test Cases (Examples provided)
# Examples of diglyceride SMILES
smiles_examples = [
    "CCCCCCCCCCCC(=O)OCC(CO)OC(=O)CCCCCCC/C=C\\CCCCCCCC",  # 1-lauroyl-2-oleoylglycerol
    "O(C(=O)CCCCCCCCC/C=C\\C/C=C\\CCCCC)C[C@@H](O)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC",  # DG(20:2n6/0:0/22:5n6)
    "O(C[C@@H](OC(CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)=O)CO)C(=O)CCCC/C=C\\C/C=C\\C/C=C\\CCCCC"  # DG(18:3(6Z,9Z,12Z)/20:4(5Z,8Z,11Z,14Z)/0:0)
]

for smiles in smiles_examples:
    result, reason = is_diglyceride(smiles)
    print(f"SMILES: {smiles}\nResult: {result}, Reason: {reason}\n")
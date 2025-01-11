"""
Classifies: CHEBI:18035 diglyceride
"""
from rdkit import Chem

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

    # Look for glycerol-like backbone (C-C-C with two -OH in any arrangement)
    glycerol_pattern = Chem.MolFromSmarts("[O][CH2][CH](O[CX3,CX2])[CH2]O[CX3,CX2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No suitable glycerol-like backbone found"
        
    # Look for exactly 2 ester groups (-O-C(=O)-) 
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0][#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Counting carbons to ensure reasonably long fatty acid chains
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10: # Assuming a minimum combined carbon count for fatty acids
        return False, "Too few carbons, likely not a diglyceride"

    return True, "Contains glycerol backbone with two fatty acid chains attached via ester bonds"

# Example Test Cases for the function
smiles_examples = [
    "CCCCCCCCCCCC(=O)OCC(CO)OC(=O)CCCCCCC/C=C\\CCCCCCCC",  # 1-lauroyl-2-oleoylglycerol
    "O(C(=O)CCCCCCCCC/C=C\\C/C=C\\CCCCC)C[C@@H](O)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC",  # DG(20:2n6/0:0/22:5n6)
    "O(C[C@@H](OC(CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)=O)CO)C(=O)CCCC/C=C\\C/C=C\\C/C=C\\CCCCC"  # DG(18:3(6Z,9Z,12Z)/20:4(5Z,8Z,11Z,14Z)/0:0)
]

for smiles in smiles_examples:
    result, reason = is_diglyceride(smiles)
    print(f"SMILES: {smiles}\nResult: {result}, Reason: {reason}\n")
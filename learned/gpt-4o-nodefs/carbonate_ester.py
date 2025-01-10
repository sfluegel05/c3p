"""
Classifies: CHEBI:46722 carbonate ester
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule is a carbonate ester based on its SMILES string.
    A carbonate ester typically contains an R-O-C(=O)-O-R' structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbonate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carbonate ester pattern (R-O-(C=O)-O-R'), allowing for various R groups
    # The SMARTS pattern [OD2] denotes a double-bonded oxygen with two connections,
    # and [OD2]C(=O)[OD2] indicates a carbon doubly bonded to one of those oxygens and singly bonded to another.
    carbonate_ester_pattern = Chem.MolFromSmarts("[$([OX2H][C]=O[OX2H1])]")
    
    if mol.HasSubstructMatch(carbonate_ester_pattern):
        return True, "Contains carbonate ester group (R-O-(C=O)-O-R')"
    else:
        return False, "No carbonate ester group found"

# Test cases with example SMILES, to validate improvement
examples = [
    "O=C1O[C@@H]2C=C[C@@H]([C@]([C@@H]2O1)(O)C#CC(=C)C)O",  # Stagonosporyne E
    "CC/C=C\\CCOC(OC)=O", # (Z)-hex-3-en-1-yl methyl carbonate
    "O=C1O[C@H]2[C@@](O)(C#CC(=C)C)[C@@H](O)CC[C@H]2O1",  # Stagonosporyne D
    "CCOC(=O)OC1=C(C(=O)N(C)C11CCN(CC1)OC)C1=C(C)C=C(Cl)C=C1C",  # spiropidion
    "CCC(C)c1cc(cc(c1OC(=O)OC(C)C)[N+]([O-])=O)[N+]([O-])=O",  # False positive example from previous false positives
    "O(CCC(C)=C)C(OCC)=O" # Ethyl 3-methylbut-3-enyl carbonate
]

for smiles in examples:
    result, reason = is_carbonate_ester(smiles)
    print(f"SMILES: {smiles} -> Is Carbonate Ester: {result}, Reason: {reason}")
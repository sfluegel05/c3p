"""
Classifies: CHEBI:46722 carbonate ester
"""
from rdkit import Chem

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule is a carbonate ester based on its SMILES string.
    A carbonate ester typically involves an R-O-C(=O)-O-R' functionality, where R and R' 
    can be any generic organic groups, potentially within rings or complex systems.

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
    
    # General carbonate ester pattern accounting for varied structure: incorporating ester (-O-C(=O)-O-)
    # The pattern might be generalized to allow for different attachment points and flexible linkage.
    carbonate_ester_pattern = Chem.MolFromSmarts("C(=O)(O*)O*")

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
    "CCCC(=C(C(=C(C(=C(CCC1=CC=C(O1)OC(=O)O)[H])[H])[H])[H])[H])[H]", # 5-(undeca-3,5,7-trien-1-yl)-2-furyl hydrogen carbonate
    "O=C1O[C@H]2C(=C(CO)[C@@H]([C@@H]([C@H]2O1)O)O)/C=C/CCCCC" # Phomoxin
]

for smiles in examples:
    result, reason = is_carbonate_ester(smiles)
    print(f"SMILES: {smiles} -> Is Carbonate Ester: {result}, Reason: {reason}")
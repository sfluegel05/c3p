"""
Classifies: CHEBI:22307 aldoxime
"""
from rdkit import Chem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime based on its SMILES string.
    An aldoxime is defined as an oxime derived from an aldehyde,
    typically characterized by the structure RCH=NOH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldoxime, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define more elaborate aldoxime patterns
    # Primary pattern for aldoxime: where =N-O must belong to RCH=NOH
    aldoxime_patterns = [
        Chem.MolFromSmarts("[#6][CH]=[N][OX1H]"), # Main definition pattern
        Chem.MolFromSmarts("[#6][CH]([#6])[CX3](=[N][OX1H])"), # Considering alkyl groups near aldoxime
    ]

    # Attempt to exclude patterns similar to historical false positives
    false_positive_patterns = [
        Chem.MolFromSmarts("[#6]=[N]/[OX2]"),  # Exclude sulfonamides and similar structures
        Chem.MolFromSmarts("[#7](=[#8])[#7]"), # Exclude N-oxide
    ]

    # Check false positive patterns first
    for pattern in false_positive_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "False positive pattern detected"

    # Match against the refined aldoxime patterns
    for pattern in aldoxime_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains aldoxime functional group"

    return False, "No aldoxime functional group found"

# Example usage with example SMILES
smiles_examples = [
    "CN1C=C(C=C1C=NO)C(=O)C2=CC=CC(=C2)Cl", # 4-(3-chlorobenzoyl)-1-methyl-pyrrole-2-carbaldehyde oxime
    "C([C@@H](/C(=N/O)/[H])C)C", # (1E,2S)-2-methylbutanal oxime
    "O=C1OC2=C(C=CC(=C2)C3=C(OC(=O)C)C(OC(=O)C)=C(C4=CC=C(OC(=O)C)C=C4)C(=C3O)O)OC1(N(O)C(=O)/C(=N/O)/[C@H](CC)C)[C@H](CC)C" # Sarcodonin, to verify exclusion
]
results = [is_aldoxime(smiles) for smiles in smiles_examples]
for (smiles, res) in zip(smiles_examples, results):
    print(f"SMILES: {smiles}, Result: {res}")
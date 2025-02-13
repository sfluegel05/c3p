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
    aldoxime_pattern = Chem.MolFromSmarts("[CH]=[N][OH]")  # Aldehyde-derived oxime pattern

    # Exclude N-oxide pattern
    false_positive_pattern = Chem.MolFromSmarts("[N]=[O]")  # Typical structure for N-oxide
    
    # Check false positive patterns first
    if mol.HasSubstructMatch(false_positive_pattern):
        return False, "False positive pattern (N-oxide) detected"

    # Match against the aldoxime pattern
    if mol.HasSubstructMatch(aldoxime_pattern):
        return True, "Contains aldoxime functional group"

    return False, "No aldoxime functional group found"

# Example usage with example SMILES
smiles_examples = [
    "CN1C=C(C=C1C=NO)C(=O)C2=CC=CC(=C2)Cl",  # 4-(3-chlorobenzoyl)-1-methyl-pyrrole-2-carbaldehyde oxime
    "C([C@@H](/C(=N/O)/[H])C)C",              # (1E,2S)-2-methylbutanal oxime
    "N(=CC1=CC=C(O1)[N+]([O-])=O)O"           # nifuroxime
]
results = [is_aldoxime(smiles) for smiles in smiles_examples]
for (smiles, res) in zip(smiles_examples, results):
    print(f"SMILES: {smiles}, Result: {res}")
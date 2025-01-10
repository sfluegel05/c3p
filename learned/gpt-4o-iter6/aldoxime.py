"""
Classifies: CHEBI:22307 aldoxime
"""
from rdkit import Chem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime based on its SMILES string.
    An aldoxime is defined as an oxime derived from an aldehyde,
    typically characterized by the structure RCHO where the CHO is converted to CH=NOH.

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
    
    # Correcting the check for aldoxime pattern
    # Search for C=N-O pattern where C is not fully saturated (aldehyde-derived)
    # Ensure we capture both possible isomerisms E/Z
    aldoxime_pattern = Chem.MolFromSmarts("[#6]=[N][OX1]")
    if mol.HasSubstructMatch(aldoxime_pattern):
        return True, "Contains aldoxime functional group"

    return False, "No aldoxime functional group found"

# Example test
smiles_examples = [
    "CN1C=C(C=C1C=NO)C(=O)C2=CC=CC(=C2)Cl", # 4-(3-chlorobenzoyl)-1-methyl-pyrrole-2-carbaldehyde oxime
    "C([C@@H](/C(=N/O)/[H])C)C" # (1E,2S)-2-methylbutanal oxime
]
results = [is_aldoxime(smiles) for smiles in smiles_examples]
for (smiles, res) in zip(smiles_examples, results):
    print(f"SMILES: {smiles}, Result: {res}")
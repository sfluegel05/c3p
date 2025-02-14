"""
Classifies: CHEBI:36141 quinone
"""
from rdkit import Chem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    Quinones are compounds with a fully conjugated cyclic dione structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for a simple archetypical quinone: conjugated dione in an aromatic ring
    quinone_pattern = Chem.MolFromSmarts("[$([cR1]=O)]~[$([cR1]=O)]")

    # Check for quinone pattern
    if mol.HasSubstructMatch(quinone_pattern):
        return True, "Contains fully conjugated cyclic dione structure typical of quinones"
    
    return False, "Does not contain the structural features typical of quinones"

# Examples
quinone_smiles = [
    "Cc1ccc2C(=O)c3ccccc3C(=O)c2c1",  # Example SMILES with a naphthoquinone structure
    "O=C1OC(C2=C(O)C(=CC3=C2C(=O)C(N)=CC3=O)C)=C(C)C=C1C"  # Salinaphthoquinone C
]

for smiles in quinone_smiles:
    result, reason = is_quinone(smiles)
    print(f"SMILES: {smiles} -> Quinone: {result}, Reason: {reason}")
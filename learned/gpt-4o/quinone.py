"""
Classifies: CHEBI:36141 quinone
"""
from rdkit import Chem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    Quinones are characterized by a fully conjugated cyclic dione structure.

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
    
    # SMARTS pattern for capturing various fully conjugated cyclic dione structures
    quinone_patterns = [
        Chem.MolFromSmarts("c1cc(=O)cc(=O)c1"),  # Simple benzoquinone-like
        Chem.MolFromSmarts("c1ccc2C(=O)cc(=O)c2c1"),  # Naphthoquinone-like
        Chem.MolFromSmarts("c1ccc2c(c1)c(=O)c3ccccc3c2=O"),  # Anthraquinone-like
    ]
    
    # Check for quinone patterns
    for pattern in quinone_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains fully conjugated cyclic dione structure typical of quinones"
    
    return False, "Does not contain the structural features typical of quinones"

# Example testing with known quinone SMILES strings
quinone_smiles = [
    "Cc1ccc2C(=O)c3ccccc3C(=O)c2c1",  # Example SMILES with naphthoquinone structure
    "O=C1OC(C2=C(O)C(=CC3=C2C(=O)C(N)=CC3=O)C)=C(C)C=C1C"  # Salinaphthoquinone C
]

for smiles in quinone_smiles:
    result, reason = is_quinone(smiles)
    print(f"SMILES: {smiles} -> Quinone: {result}, Reason: {reason}")
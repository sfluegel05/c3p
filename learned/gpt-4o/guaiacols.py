"""
Classifies: CHEBI:134251 guaiacols
"""
from rdkit import Chem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    A guaiacol should have a phenol group with an additional methoxy substituent at the ortho-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a guaiacol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More flexible guaiacol SMARTS pattern to match ortho methoxy and phenol
    # Coving any substitution around the benzene, checking only ortho methoxy and phenol
    flexible_guaiacol_pattern = Chem.MolFromSmarts("Oc1cc(ccc1)CO")

    # Check for the guaiacol pattern
    if mol.HasSubstructMatch(flexible_guaiacol_pattern):
        return True, "Molecule contains a phenol group with an ortho methoxy substituent"

    return False, "Molecule does not match the guaiacol structure with phenol and ortho methoxy"

# Example usage
smiles_examples = [
    "COc1ccccc1O",  # guaiacol
    "C1=CC=C(C=C1O)C=C(O)COC"  # Example of similar structures for clarity
]

for smiles in smiles_examples:
    result, reason = is_guaiacols(smiles)
    print(f"SMILES: {smiles} => Guaiacol: {result}. Reason: {reason}")
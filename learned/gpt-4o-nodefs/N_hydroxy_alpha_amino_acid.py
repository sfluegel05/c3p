"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
from rdkit import Chem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-hydroxy-alpha-amino-acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-hydroxy-alpha-amino-acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the key structural patterns for N-hydroxy-alpha-amino-acid
    
    # This SMARTS pattern matches an alpha amino acid backbone with a hydroxyamino modification:
    nhydroxy_amino_pattern = Chem.MolFromSmarts("N(O)[C@@H]([*!#0])[CX3](=O)[OX1H]")
    if not mol.HasSubstructMatch(nhydroxy_amino_pattern):
        return False, "Does not contain the N-hydroxy alpha-amino acid substructure"
    
    return True, "Contains N-hydroxy alpha-amino acid substructure"

# Test examples: These should return True
example_smiles = [
    r"O=C(O)[C@@H](N(O)O)CCCCCCCSC",
    r"C(\N)(=N/O)/NCCC[C@H](N)C(=O)O",
    r"NC(CCCNC(=N)NO)C(O)=O",
    # Add more test examples as needed
]

for smile in example_smiles:
    result, reason = is_N_hydroxy_alpha_amino_acid(smile)
    print(f"SMILES: {smile}, Is N-hydroxy-alpha-amino-acid: {result}, Reason: {reason}")
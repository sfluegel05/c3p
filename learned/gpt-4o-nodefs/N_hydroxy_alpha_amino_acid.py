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
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the key structural patterns for N-hydroxy-alpha-amino-acid
    # The pattern captures an amino group with an attached hydroxy group (N-O-H)
    # and an alpha carbon with a carboxyl group (C(=O)O)
    nhydroxy_alpha_amino_pattern = Chem.MolFromSmarts("[CX4][CX3](N(O))(*)[CX3](=O)[OX1H]")
    
    if not mol.HasSubstructMatch(nhydroxy_alpha_amino_pattern):
        return False, "Does not contain N-hydroxy-alpha-amino-acid substructure"
    
    return True, "Contains N-hydroxy-alpha-amino-acid substructure"

# Test examples: These should return True
example_smiles = [
    "O=C(O)[C@@H](N(O)O)CCCCCCCSC",
    "C(\\N)(=N/O)/NCCC[C@H](N)C(=O)O",
    "NC(CCCNC(=N)NO)C(O)=O",
    "O=C(O)[C@@H](N(O)O)CCCCCCCCSC",
    "N1([C@@H](CCCC1)C(=O)O)O",
    "ON(O)[C@@H](Cc1ccccc1)C(O)=O",
]

for smile in example_smiles:
    result, reason = is_N_hydroxy_alpha_amino_acid(smile)
    print(f"SMILES: {smile}, Is N-hydroxy-alpha-amino-acid: {result}, Reason: {reason}")
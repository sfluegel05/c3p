"""
Classifies: CHEBI:63534 monoamine
"""
"""
Classifies: MONOAMINE
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine is an aralylamino compound which contains one amino group connected to 
    an aromatic ring by a two-carbon chain. Monoamines are derived from aromatic amino acids 
    like phenylalanine, tyrosine, tryptophan, and the thyroid hormones by the action of 
    aromatic amino acid decarboxylase enzymes.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for amino group connected via a two-carbon chain to an aromatic ring
    monoamine_pattern = Chem.MolFromSmarts("[N;!R;X3]-[C;!R;X4]-[C;!R;X4]-[a]")
    
    # Find matches of the pattern in the molecule
    matches = mol.GetSubstructMatches(monoamine_pattern)
    
    if len(matches) == 1:
        return True, "Contains one amino group connected via two-carbon chain to an aromatic ring"
    elif len(matches) == 0:
        return False, "No amino group connected via two-carbon chain to an aromatic ring found"
    else:
        return False, f"Found {len(matches)} amino groups connected via two-carbon chains to aromatic rings, expected exactly one"
"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
from rdkit import Chem

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone based on its SMILES string.
    A hydroxynaphthoquinone is a naphthoquinone moiety with at least one hydroxy group substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxynaphthoquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Refined SMARTS pattern for naphthoquinone
    naphthoquinone_patterns = [
        Chem.MolFromSmarts('c1ccccc2C(=O)C=CC2=O'), # Base pattern for naphthoquinone
        Chem.MolFromSmarts('C1=CC=C2C(=O)C=CC2=O')  # Alternate structure
    ]
    hydroxy_pattern = Chem.MolFromSmarts('[OH]')

    # Check for naphthoquinone structure
    for pattern in naphthoquinone_patterns:
        if mol.HasSubstructMatch(pattern):
            naphtho_matches = mol.GetSubstructMatches(pattern)
            break
    else:
        return False, "No naphthoquinone core found"

    # Check for hydroxy group in relation to the naphthoquinone
    for naphtho_match in naphtho_matches:
        naphtho_atoms = set(naphtho_match)
        for atom_idx in naphtho_match:
            atom_obj = mol.GetAtomWithIdx(atom_idx)
            if atom_obj.GetSymbol() == 'C':  # Look for C atoms to attach OH
                neighbors = atom_obj.GetNeighbors()
                for neighbor in neighbors:
                    if neighbor.GetSymbol() == 'O':
                        for o_neighbor in neighbor.GetNeighbors():
                            if o_neighbor.GetIdx() not in naphtho_atoms:
                                return True, "Contains naphthoquinone core with hydroxy group substitution"

    return False, "No correctly attached hydroxy group found"
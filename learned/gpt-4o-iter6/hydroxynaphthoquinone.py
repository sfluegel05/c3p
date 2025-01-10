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

    # Define SMARTS patterns - more flexible pattern for naphthoquinone
    naphthoquinone_pattern = Chem.MolFromSmarts('c1ccc2c(c1)C(=O)C=CC2=O')
    hydroxy_pattern = Chem.MolFromSmarts('[OH]')

    # Check for naphthoquinone structure
    if not mol.HasSubstructMatch(naphthoquinone_pattern):
        return False, "No naphthoquinone core found"

    # Check for at least one hydroxy group attached to the aromatic system
    naphtho_matches = mol.GetSubstructMatches(naphthoquinone_pattern)
    for naphtho_match in naphtho_matches:
        naphtho_atoms = set(naphtho_match)
        for atom in naphtho_match:
            atom_obj = mol.GetAtomWithIdx(atom)
            if atom_obj.GetSymbol() == 'C':  # Look for C atoms to attach OH
                neighbors = atom_obj.GetNeighbors()
                for neighbor in neighbors:
                    if neighbor.GetSymbol() == 'O' and neighbor.GetIdx() not in naphtho_atoms:
                        return True, "Contains naphthoquinone core with hydroxy group substitution"

    return False, "Hydroxy group not attached correctly to naphthoquinone"
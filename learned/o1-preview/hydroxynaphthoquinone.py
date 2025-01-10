"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
"""
Classifies: CHEBI:51793 hydroxynaphthoquinone
"""
from rdkit import Chem

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone based on its SMILES string.
    A hydroxynaphthoquinone is any naphthoquinone where the naphthoquinone moiety is substituted by at least one hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxynaphthoquinone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for naphthoquinone cores (1,4- and 1,2-naphthoquinone)
    naphthoquinone_smarts_list = [
        'C12=C(C=CC=C1)C(=O)C=CC2=O',   # 1,4-naphthoquinone core
        'C1=CC=C2C(=O)C=CC(=O)C2=C1'    # 1,2-naphthoquinone core
    ]
    naphthoquinone_mols = [Chem.MolFromSmarts(smarts) for smarts in naphthoquinone_smarts_list]

    # Check for presence of naphthoquinone moiety
    has_naphthoquinone = False
    for naphthoquinone in naphthoquinone_mols:
        if mol.HasSubstructMatch(naphthoquinone):
            has_naphthoquinone = True
            naphthoquinone_matches = mol.GetSubstructMatches(naphthoquinone)
            break

    if not has_naphthoquinone:
        return False, "No naphthoquinone moiety found"

    # Check for hydroxy substitution on the naphthoquinone ring
    for match in naphthoquinone_matches:
        naphthoquinone_atoms = set(match)
        has_hydroxy = False

        for idx in naphthoquinone_atoms:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6 and atom.IsInRing():
                for neighbor in atom.GetNeighbors():
                    # Look for hydroxy group attached to ring carbon
                    if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                        has_hydroxy = True
                        break
                if has_hydroxy:
                    break

        if has_hydroxy:
            return True, "Contains naphthoquinone moiety substituted with at least one hydroxy group"

    return False, "No hydroxy substitution on naphthoquinone moiety found"
           
__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:51793',
        'name': 'hydroxynaphthoquinone',
        'definition': 'Any naphthoquinone in which the naphthoquinone moiety is substituted by at least one hydroxy group.',
        'parents': ['CHEBI:36124']  # CHEBI:36124 is naphthoquinone
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}
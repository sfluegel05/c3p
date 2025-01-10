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

    # Define SMARTS pattern for naphthoquinone core (both 1,4- and 1,2-naphthoquinone)
    naphthoquinone_smarts = Chem.MolFromSmarts('C1=CC=CC=C1C(=O)C=CC=O')  # Simplified pattern
    if not mol.HasSubstructMatch(naphthoquinone_smarts):
        return False, "No naphthoquinone core found"

    # Get the matches of the naphthoquinone core
    matches = mol.GetSubstructMatches(naphthoquinone_smarts)
    if not matches:
        return False, "No naphthoquinone core found"

    # Define SMARTS pattern for hydroxy group attached to an aromatic ring carbon
    hydroxy_smarts = Chem.MolFromSmarts('[OX1H]-[#6]')  # Hydroxy group attached to carbon

    # Check if any atom in the naphthoquinone core is substituted with a hydroxy group
    for match in matches:
        naphthoquinone_atoms = set(match)
        # Check all atoms for hydroxy substitution
        for atom_idx in naphthoquinone_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Skip non-carbon atoms
            if atom.GetAtomicNum() != 6:
                continue
            # Check neighbors for hydroxy group
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalDegree() == 1:
                    bond = mol.GetBondBetweenAtoms(atom_idx, neighbor.GetIdx())
                    if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                        return True, "Contains naphthoquinone moiety substituted with at least one hydroxy group"

    return False, "Naphthoquinone core found but no hydroxy substitution on the ring carbons"

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
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}
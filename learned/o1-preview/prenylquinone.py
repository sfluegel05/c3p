"""
Classifies: CHEBI:26255 prenylquinone
"""
"""
Classifies: prenylquinone
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    A prenylquinone is a quinone substituted by a polyprenyl-derived side-chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenylquinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for quinone ring (1,4-benzoquinone and 1,4-naphthoquinone)
    quinone_smarts = ['O=C1C=CC=CC1=O', 'O=C1C=CC2=CC=CC=C12']  # benzoquinone and naphthoquinone
    quinone_patterns = [Chem.MolFromSmarts(s) for s in quinone_smarts]

    # Search for quinone ring
    quinone_matches = []
    for qp in quinone_patterns:
        matches = mol.GetSubstructMatches(qp)
        if matches:
            quinone_matches.extend(matches)
    if not quinone_matches:
        return False, "No quinone ring found"

    # Define SMARTS pattern for prenyl unit
    prenyl_smarts = '[CH2]=C([CH3])[CH2]'  # Prenyl group
    prenyl_pattern = Chem.MolFromSmarts(prenyl_smarts)

    # Identify atoms in quinone ring(s)
    quinone_ring_atoms = set()
    for match in quinone_matches:
        quinone_ring_atoms.update(match)

    # Check for polyprenyl-derived side chain attached to quinone ring
    side_chain_found = False
    for atom_idx in quinone_ring_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for bond in atom.GetBonds():
            neighbor = bond.GetOtherAtom(atom)
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in quinone_ring_atoms:
                # Found an atom outside the quinone ring
                # Check if this leads to a polyprenyl side chain
                side_chain_atoms = set()
                atoms_to_visit = [neighbor_idx]
                while atoms_to_visit:
                    current_idx = atoms_to_visit.pop()
                    if current_idx not in side_chain_atoms:
                        side_chain_atoms.add(current_idx)
                        current_atom = mol.GetAtomWithIdx(current_idx)
                        for nbr in current_atom.GetNeighbors():
                            nbr_idx = nbr.GetIdx()
                            if nbr_idx not in side_chain_atoms:
                                atoms_to_visit.append(nbr_idx)
                # Create sub-molecule of the side chain
                side_chain_mol = Chem.PathToSubmol(mol, list(side_chain_atoms))
                # Search for prenyl repeats
                num_prenyl_units = len(side_chain_mol.GetSubstructMatches(prenyl_pattern))
                if num_prenyl_units >= 2:
                    side_chain_found = True
                    break
        if side_chain_found:
            break

    if not side_chain_found:
        return False, "No polyprenyl-derived side chain attached to quinone ring"

    return True, "Contains quinone ring with polyprenyl-derived side chain"

__metadata__ = {
    'chemical_class': {
        'id': '',
        'name': 'prenylquinone',
        'definition': 'A quinone substituted by a polyprenyl-derived side-chain. Prenylquinones occur in all living cells. Due to their amphiphilic character, they are mainly located in biological membranes where they function as electron and proton carriers in the photosynthetic and respiratory electron transport chains. Some prenylquinones also perform more specialised roles such as antioxidants and enzyme cofactors. Prenylquinones are classified according to ring structure: the main classes are menaquinones, phylloquinones, ubiquinones and plastoquinones.',
        'parents': []
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
    'stdout': None,
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}
"""
Classifies: CHEBI:84403 phytoceramides
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phytoceramides(smiles: str):
    """
    Determines if a molecule is a phytoceramide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phytoceramide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a sphingoid base hydroxylated at position 4
    hydroxylated_sphingoid = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and len(atom.GetNeighbors()) == 1:
            neighbor = atom.GetNeighbors()[0]
            if neighbor.GetSymbol() == 'C' and len(neighbor.GetNeighbors()) == 4:
                hydroxylated_sphingoid = True
                break

    if not hydroxylated_sphingoid:
        return False, "No hydroxylated sphingoid base at position 4 found"

    # Check for the presence of long carbon chains (14 to 20 carbons) in the backbone
    carbon_chains = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and len(atom.GetNeighbors()) == 2:
            chain_length = 1
            current_atom = atom
            while len(current_atom.GetNeighbors()) == 2:
                next_atom = [neighbor for neighbor in current_atom.GetNeighbors() if neighbor.GetSymbol() == 'C' and neighbor != atom][0]
                chain_length += 1
                current_atom = next_atom
            if 14 <= chain_length <= 20:
                carbon_chains.append(chain_length)

    if not carbon_chains:
        return False, "No long carbon chains (14 to 20 carbons) found in the backbone"

    # Check for the presence of N-linked fatty acids with chain lengths from 14 to 26 carbons
    n_linked_fatty_acids = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C' and len(neighbor.GetNeighbors()) == 2:
                    chain_length = 1
                    current_atom = neighbor
                    while len(current_atom.GetNeighbors()) == 2:
                        next_atom = [n for n in current_atom.GetNeighbors() if n.GetSymbol() == 'C' and n != atom][0]
                        chain_length += 1
                        current_atom = next_atom
                    if 14 <= chain_length <= 26:
                        n_linked_fatty_acids.append(chain_length)

    if not n_linked_fatty_acids:
        return False, "No N-linked fatty acids with chain lengths from 14 to 26 carbons found"

    # Check for the presence of a hydroxyl group on carbon 2
    hydroxyl_on_c2 = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and len(atom.GetNeighbors()) == 4:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O' and len(neighbor.GetNeighbors()) == 1:
                    hydroxyl_on_c2 = True
                    break

    if not hydroxyl_on_c2:
        return False, "No hydroxyl group on carbon 2 found"

    return True, "Molecule is a phytoceramide"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:84403',
                          'name': 'phytoceramides',
                          'definition': 'A phytoceramide is a ceramide where '
                                        'the sphingoid base is hydroxylated at '
                                        'position 4. These backbone bases can '
                                        'be 14 to 20 carbons long, either '
                                        'straight chain or branched. The '
                                        'N-linked fatty acids are typically '
                                        'saturated or monounsaturated with '
                                        'chain lengths from 14 to 26 carbon '
                                        'atoms; the presence of a hydroxyl '
                                        'group on carbon 2 is fairly common. '
                                        'In the illustrated generalised '
                                        'structure, R1 = OH, OX (where X = '
                                        'acyl, glycosyl, phosphate, '
                                        'phosphonate, etc.), or H.',
                          'parents': ['CHEBI:17761']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 17,
    'num_false_negatives': 17,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}
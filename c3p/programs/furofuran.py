"""
Classifies: CHEBI:47790 furofuran
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_furofuran(smiles: str):
    """
    Determines if a molecule is a furofuran (Organic heterobicyclic compounds containing two furan rings ortho-fused to each other).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furofuran, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for at least two 5-membered rings
    five_membered_rings = [ring for ring in rings.AtomRings() if len(ring) == 5]
    if len(five_membered_rings) < 2:
        return False, "Less than two 5-membered rings found"

    # Find all furan rings (5-membered rings containing one oxygen and four carbons)
    furan_rings = []
    for ring in five_membered_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if sum(1 for atom in atoms if atom.GetSymbol() == 'O') == 1 and sum(1 for atom in atoms if atom.GetSymbol() == 'C') == 4:
            furan_rings.append(ring)

    if len(furan_rings) < 2:
        return False, "Less than two furan rings found"

    # Check if the furan rings are ortho-fused
    def is_ortho_fused(ring1, ring2):
        shared_atoms = set(ring1) & set(ring2)
        return len(shared_atoms) == 2

    for i in range(len(furan_rings)):
        for j in range(i + 1, len(furan_rings)):
            if is_ortho_fused(furan_rings[i], furan_rings[j]):
                return True, "Molecule is a furofuran"

    return False, "No ortho-fused furan rings found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:47790',
                          'name': 'furofuran',
                          'definition': 'Organic heterobicyclic compounds '
                                        'containing a two furan rings '
                                        'ortho-fused to each other.',
                          'parents': ['CHEBI:27171', 'CHEBI:38104']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 15,
    'num_false_positives': 0,
    'num_true_negatives': 15,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}
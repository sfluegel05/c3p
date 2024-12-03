"""
Classifies: CHEBI:37581 gamma-lactone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma-lactone (a lactone having a five-membered lactone ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a gamma-lactone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for at least one 5-membered ring
    if not any(len(ring) == 5 for ring in rings.AtomRings()):
        return False, "No 5-membered rings found"

    # Find all 5-membered rings
    five_membered_rings = [ring for ring in rings.AtomRings() if len(ring) == 5]

    # Check if any of the 5-membered rings is a lactone
    for ring in five_membered_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        # Check for presence of one oxygen and one carbonyl group (C=O) in the ring
        oxygens = [atom for atom in atoms if atom.GetSymbol() == 'O']
        carbonyl_carbons = [atom for atom in atoms if atom.GetSymbol() == 'C' and any(bond.GetBondTypeAsDouble() == 2.0 and mol.GetAtomWithIdx(bond.GetOtherAtomIdx(atom.GetIdx())).GetSymbol() == 'O' for bond in atom.GetBonds())]

        if len(oxygens) == 1 and len(carbonyl_carbons) == 1:
            return True, "Gamma-lactone found"

    return False, "No gamma-lactone ring found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37581',
                          'name': 'gamma-lactone',
                          'definition': 'A lactone having a five-membered '
                                        'lactone ring.',
                          'parents': ['CHEBI:25000']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[13:55:26] SMILES Parse Error: unclosed ring for input: '
             "'[H][C@]12OC(=O)C(=C)[C@]1([H])[C@H](O)C\\C(C)=C\\CC\\C(C)=C\x02'\n"
             '[13:55:27] SMILES Parse Error: unclosed ring for input: '
             "'OC(=O)\\C=C1OC(=O)C(Cl)=C\x01'\n",
    'stdout': '',
    'num_true_positives': 119,
    'num_false_positives': 4,
    'num_true_negatives': 16,
    'num_false_negatives': 3,
    'precision': 0.967479674796748,
    'recall': 0.9754098360655737,
    'f1': 0.9714285714285714,
    'accuracy': None}
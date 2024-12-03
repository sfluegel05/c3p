"""
Classifies: CHEBI:18154 polysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polysaccharide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule has more than ten monosaccharide residues
    ring_info = mol.GetRingInfo()
    ring_atoms = ring_info.AtomRings()
    
    monosaccharide_count = 0
    for ring in ring_atoms:
        if len(ring) == 6:  # Assuming a monosaccharide has a 6-membered ring
            monosaccharide_count += 1

    if monosaccharide_count > 10:
        return True, f"Polysaccharide with {monosaccharide_count} monosaccharide residues"
    else:
        return False, f"Only {monosaccharide_count} monosaccharide residues found, less than 11"

# Example usage
smiles = "O([C@H]1O[C@@H]([C@@H](O[C@@H]2O[C@@H]([C@@H](O[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO[C@]4(O[C@H]([C@H](NC(=O)C)[C@@H](O)C4)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O)[C@@H]1O[C@@H]5O[C@@H]([C@@H](O[C@@H]6O[C@@H]([C@H](O)[C@H](O)[C@H]6O)CO[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]5NC(=O)C)CO)CO)[C@H]8[C@H](O)[C@H](O[C@@H](O[C@H]9[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]9CO)O[C@H]%10[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]%10CO[C@@H]%11O[C@H]([C@@H](O)[C@@H](O)[C@@H]%11O)C)O)[C@H]8O)CO[C@H]%12O[C@@H]([C@@H](O)[C@H](O)[C@@H]%12O[C@@H]%13O[C@@H]([C@@H](O[C@@H]%14O[C@@H]([C@H](O)[C@H](O)[C@H]%14O)CO[C@]%15(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%15)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]%13NC(=O)C)CO)CO[C@@H]%16O[C@@H]([C@@H](O[C@@H]%17O[C@@H]([C@H](O)[C@H](O)[C@H]%17O)CO[C@]%18(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%18)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]%16NC(=O)C)CO"
print(is_polysaccharide(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:18154',
                          'name': 'polysaccharide',
                          'definition': 'A biomacromolecule consisting of '
                                        'large numbers of monosaccharide '
                                        'residues linked glycosidically. This '
                                        'term is commonly used only for those '
                                        'containing more than ten '
                                        'monosaccharide residues.',
                          'parents': [   'CHEBI:16646',
                                         'CHEBI:167559',
                                         'CHEBI:33694']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(True, 'Polysaccharide with 18 monosaccharide residues')\n",
    'num_true_positives': 164,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 1,
    'precision': 1.0,
    'recall': 0.9939393939393939,
    'f1': 0.9969604863221885,
    'accuracy': None}
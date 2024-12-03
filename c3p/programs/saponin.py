"""
Classifies: CHEBI:26605 saponin
"""
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_saponin(smiles: str):
    """
    Determines if a molecule is a saponin (a glycoside combined with a triterpenoid or steroid derivative).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saponin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycoside moiety (presence of sugar units)
    def has_glycoside(mol):
        sugars = ['C1OC(CO)C(O)C(O)C1O', 'C1OC(CO)C(O)C(O)C1O', 'C1OC(CO)C(O)C(O)C1O']  # Glucose, Galactose, Rhamnose, etc.
        for sugar in sugars:
            sugar_mol = Chem.MolFromSmiles(sugar)
            if mol.HasSubstructMatch(sugar_mol):
                return True
        return False

    if not has_glycoside(mol):
        return False, "No glycoside moiety found"

    # Check for triterpenoid or steroid backbone
    triterpenoid_steroid_motifs = [
        'C1CCC2C(C1)CCC3C2CCC4C3(CCC5C4CCC6C5(CCC(C6)O)C)C',  # Example steroid backbone
        'C1CCC2C(C1)CCC3C2CCC4C3(CCC5C4CCC6C5(CCC(C6)O)C)C'   # Example triterpenoid backbone
    ]
    has_backbone = False
    for motif in triterpenoid_steroid_motifs:
        motif_mol = Chem.MolFromSmiles(motif)
        if mol.HasSubstructMatch(motif_mol):
            has_backbone = True
            break

    if not has_backbone:
        return False, "No triterpenoid or steroid backbone found"

    return True, "Molecule is a saponin"

# Example usage:
# smiles = "O([C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@@H]([C@H](O)[C@H](O)[C@H]2O[C@@H]3O[C@H]([C@H](O)[C@@H](O)[C@H]3O)C)CO)C(O)=O)[C@@H]4[C@]([C@]5([C@](CC4)([C@@]6([C@@](CC5)([C@]7(C(=CC6)[C@]8([C@@](CC7)([C@@H]([C@@H](C(C8)(C)C)O)O)C)[H])C)C)[H])C)[H])(C)CO"
# print(is_saponin(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26605',
                          'name': 'saponin',
                          'definition': 'A glycoside that is a compound '
                                        'containing one or more hydrophilic '
                                        'glycoside moieties combined with a '
                                        'lipophilic triterpenoid or steroid '
                                        'derivative. Found in particular '
                                        'abundance in plant species.',
                          'parents': ['CHEBI:24400']},
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
    'num_true_negatives': 20,
    'num_false_negatives': 82,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}
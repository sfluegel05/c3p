"""
Classifies: CHEBI:26739 sphingolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sphingolipid(smiles: str):
    """
    Determines if a molecule is a sphingolipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingolipid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify sphingoid base backbone (sphingosine-like structure)
    sphingoid_base = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() == 'C']
            if len(neighbors) == 2:
                if any(n.GetSymbol() == 'O' for n in neighbors[0].GetNeighbors()) and \
                   any(n.GetSymbol() == 'O' for n in neighbors[1].GetNeighbors()):
                    sphingoid_base = True
                    break

    if not sphingoid_base:
        return False, "No sphingoid base backbone found"

    # Check for long-chain fatty acid amide linkage
    fatty_acid_amide = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and \
           bond.GetBeginAtom().GetSymbol() == 'N' and \
           bond.GetEndAtom().GetSymbol() == 'C':
            carbon = bond.GetEndAtom() if bond.GetBeginAtom().GetSymbol() == 'N' else bond.GetBeginAtom()
            if carbon.GetTotalNumHs() == 2 and carbon.GetDegree() == 3:
                fatty_acid_amide = True
                break

    if not fatty_acid_amide:
        return False, "No long-chain fatty acid amide linkage found"

    return True, "Molecule is a sphingolipid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26739',
                          'name': 'sphingolipid',
                          'definition': 'Sphingolipids are a complex family of '
                                        'compounds that share a common '
                                        'structural feature, a sphingoid base '
                                        'backbone.',
                          'parents': [   'CHEBI:18059',
                                         'CHEBI:35352',
                                         'CHEBI:36963']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[15:23:43] Explicit valence for atom # 1 Cl, 4, is greater than '
             'permitted\n',
    'stdout': '',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 229,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}
"""
Classifies: CHEBI:35627 beta-lactam
"""
from rdkit import Chem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    rings = mol.GetRingInfo()

    # Check for 4-membered rings
    four_membered_rings = [ring for ring in rings.AtomRings() if len(ring) == 4]
    if not four_membered_rings:
        return False, "No 4-membered rings found"

    for ring in four_membered_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        # Check if the ring contains an amide bond (C=O and N)
        contains_carbonyl = any(atom.GetSymbol() == 'C' and any(bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and mol.GetAtomWithIdx(bond.GetOtherAtomIdx(atom.GetIdx())).GetSymbol() == 'O' for bond in atom.GetBonds()) for atom in atoms)
        contains_amide_nitrogen = any(atom.GetSymbol() == 'N' and any(bond.GetBondType() == Chem.rdchem.BondType.SINGLE and mol.GetAtomWithIdx(bond.GetOtherAtomIdx(atom.GetIdx())).GetSymbol() == 'C' for bond in atom.GetBonds()) for atom in atoms)
        
        if contains_carbonyl and contains_amide_nitrogen:
            return True, "Contains a 4-membered ring with an amide bond"

    return False, "No 4-membered rings with an amide bond found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35627',
                          'name': 'beta-lactam',
                          'definition': 'A lactam in which the amide bond is '
                                        'contained within a four-membered '
                                        'ring, which includes the amide '
                                        'nitrogen and the carbonyl carbon.',
                          'parents': ['CHEBI:24995']},
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
    'num_true_positives': 33,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 1,
    'precision': 1.0,
    'recall': 0.9705882352941176,
    'f1': 0.9850746268656716,
    'accuracy': None}
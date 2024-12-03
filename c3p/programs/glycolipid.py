"""
Classifies: CHEBI:33563 glycolipid
"""
from rdkit import Chem

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycolipid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of glycosidic linkage (C-O-C bond)
    glycosidic_linkage_found = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if (begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 8) or (begin_atom.GetAtomicNum() == 8 and end_atom.GetAtomicNum() == 6):
                for neighbor in end_atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != begin_atom.GetIdx():
                        glycosidic_linkage_found = True
                        break
                if glycosidic_linkage_found:
                    break

    if not glycosidic_linkage_found:
        return False, "No glycosidic linkage found"

    # Check for the presence of carbohydrate part (ring structure with oxygen and multiple hydroxyl groups)
    carbohydrate_found = False
    for ring in mol.GetRingInfo().AtomRings():
        oxygens = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        hydroxyls = sum(1 for idx in ring if any(neighbor.GetAtomicNum() == 8 for neighbor in mol.GetAtomWithIdx(idx).GetNeighbors()))
        if oxygens >= 1 and hydroxyls >= 2:
            carbohydrate_found = True
            break

    if not carbohydrate_found:
        return False, "No carbohydrate part found"

    # Check for the presence of 1,2-di-O-acylglycerol (glycerol backbone with two acyl groups)
    diacylglycerol_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetDegree() == 3:
            oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
            if len(oxygens) == 2:
                acyl_groups = sum(1 for o in oxygens if any(n.GetAtomicNum() == 6 and n.GetDegree() == 3 for n in o.GetNeighbors()))
                if acyl_groups == 2:
                    diacylglycerol_found = True
                    break

    if not diacylglycerol_found:
        return False, "No 1,2-di-O-acylglycerol found"

    return True, "Molecule is a glycolipid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33563',
                          'name': 'glycolipid',
                          'definition': 'Any member of class of '
                                        '1,2-di-O-acylglycerols joined at '
                                        'oxygen 3 by a glycosidic linkage to a '
                                        'carbohydrate part (usually a mono-, '
                                        'di- or tri-saccharide). Some '
                                        'substances classified as bacterial '
                                        'glycolipids have the sugar part '
                                        'acylated by one or more fatty acids '
                                        'and the glycerol part may be absent.',
                          'parents': ['CHEBI:35740']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 161,
    'num_false_positives': 3,
    'num_true_negatives': 17,
    'num_false_negatives': 0,
    'precision': 0.9817073170731707,
    'recall': 1.0,
    'f1': 0.9907692307692307,
    'accuracy': None}
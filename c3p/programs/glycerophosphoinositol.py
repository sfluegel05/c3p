"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphoinositol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for inositol ring (a six-membered ring with five hydroxyl groups and one ether)
    inositol_found = False
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if all(atom.GetSymbol() == 'C' for atom in atoms):
                hydroxyl_count = sum(1 for atom in atoms if any(neighbor.GetSymbol() == 'O' and neighbor.GetTotalNumHs() == 1 for neighbor in atom.GetNeighbors()))
                if hydroxyl_count == 5:
                    inositol_found = True
                    break

    if not inositol_found:
        return False, "No inositol ring found"

    # Check for glycerol backbone
    glycerol_backbone_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetTotalNumHs() == 2:
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 2 and all(neighbor.GetSymbol() == 'C' and neighbor.GetTotalNumHs() == 1 for neighbor in neighbors):
                glycerol_backbone_found = True
                break

    if not glycerol_backbone_found:
        return False, "No glycerol backbone found"

    # Check for phosphate group attached to glycerol backbone
    phosphate_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'P' and atom.GetTotalNumHs() == 0:
            neighbors = atom.GetNeighbors()
            if any(neighbor.GetSymbol() == 'O' and neighbor.GetTotalNumHs() == 0 for neighbor in neighbors):
                phosphate_found = True
                break

    if not phosphate_found:
        return False, "No phosphate group found"

    # Check for ester linkage to inositol
    ester_linkage_found = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if (begin_atom.GetSymbol() == 'C' and end_atom.GetSymbol() == 'O' and end_atom.GetTotalNumHs() == 0) or \
               (end_atom.GetSymbol() == 'C' and begin_atom.GetSymbol() == 'O' and begin_atom.GetTotalNumHs() == 0):
                ester_linkage_found = True
                break

    if not ester_linkage_found:
        return False, "No ester linkage to inositol found"

    return True, "Molecule is a glycerophosphoinositol"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36315',
                          'name': 'glycerophosphoinositol',
                          'definition': 'Any glycerophospholipid having the '
                                        'polar alcohol inositol esterified to '
                                        'the phosphate group at the sn-3 '
                                        'position of the glycerol backbone.',
                          'parents': ['CHEBI:37739']},
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
    'num_true_positives': 32,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 32,
    'precision': 1.0,
    'recall': 0.5,
    'f1': 0.6666666666666666,
    'accuracy': None}
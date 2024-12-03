"""
Classifies: CHEBI:60334 peptide anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_peptide_anion(smiles: str):
    """
    Determines if a molecule is a peptide anion (An anion formed by deprotonation of at least one peptide carboxy group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a peptide anion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of peptide bonds (amide bonds)
    peptide_bonds = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if (begin_atom.GetSymbol() == 'C' and end_atom.GetSymbol() == 'N') or (begin_atom.GetSymbol() == 'N' and end_atom.GetSymbol() == 'C'):
                if begin_atom.GetDegree() > 1 and end_atom.GetDegree() > 1:
                    peptide_bonds = True
                    break

    if not peptide_bonds:
        return False, "No peptide bonds found"

    # Check for deprotonated carboxy groups
    carboxy_groups = 0
    deprotonated_carboxy_groups = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetDegree() == 3:
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            if neighbors.count('O') == 2 and neighbors.count('C') == 1:
                carboxy_groups += 1
                if any(n.GetFormalCharge() == -1 for n in atom.GetNeighbors() if n.GetSymbol() == 'O'):
                    deprotonated_carboxy_groups += 1

    if deprotonated_carboxy_groups == 0:
        return False, "No deprotonated carboxy groups found"

    return True, f"Peptide anion with {deprotonated_carboxy_groups} deprotonated carboxy group(s)"

# Example usage:
# is_peptide_anion("[NH3+][C@H](C([O-])=O)CCC(=O)NCC(=O)[O-]")
# is_peptide_anion("C(CC[C@@H](C(N1[C@@H](CCC1)C(=O)N[C@@H]([C@H](CC)C)C([O-])=O)=O)NC([C@H](C(C)C)NC([C@H](CC2=CC=CC=C2)NC(=O)[C@H]3N(CCC3)C([C@H](CC4=CC=C(C=C4)O)[NH3+])=O)=O)=O)(=O)[O-]")



__metadata__ = {   'chemical_class': {   'id': 'CHEBI:60334',
                          'name': 'peptide anion',
                          'definition': 'An anion formed by deprotonation of '
                                        'at least one peptide carboxy group.',
                          'parents': ['CHEBI:25696']},
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
    'num_true_positives': 10,
    'num_false_positives': 0,
    'num_true_negatives': 10,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}
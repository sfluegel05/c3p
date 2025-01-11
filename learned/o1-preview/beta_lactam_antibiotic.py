"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
"""
Classifies: beta-lactam antibiotic
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.
    A beta-lactam antibiotic is an organonitrogen heterocyclic antibiotic that contains a beta-lactam ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam antibiotic, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all 4-membered rings
    ssr = Chem.GetSymmSSSR(mol)
    beta_lactam_found = False
    for ring in ssr:
        if len(ring) == 4:
            atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
            atom_types = [atom.GetAtomicNum() for atom in atoms_in_ring]
            # Check if ring contains 1 nitrogen and 3 carbons
            n_count = atom_types.count(7)  # Atomic number 7 for nitrogen
            c_count = atom_types.count(6)  # Atomic number 6 for carbon
            if n_count == 1 and c_count == 3:
                # Check for carbonyl group attached to a ring carbon
                carbonyl_found = False
                for idx in ring:
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.GetAtomicNum() == 6:  # Carbon atom
                        for neighbor in atom.GetNeighbors():
                            if neighbor.GetAtomicNum() == 8:  # Oxygen atom
                                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                    carbonyl_found = True
                                    break
                        if carbonyl_found:
                            break
                if carbonyl_found:
                    beta_lactam_found = True
                    break

    if not beta_lactam_found:
        return False, "No beta-lactam ring found"

    return True, "Contains a beta-lactam ring"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'beta-lactam antibiotic',
        'definition': 'An organonitrogen heterocyclic antibiotic that contains a beta-lactam ring.',
        'parents': []
    },
    'config': {},
    'message': None,
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}
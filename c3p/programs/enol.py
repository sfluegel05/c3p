"""
Classifies: CHEBI:33823 enol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_enol(smiles: str):
    """
    Determines if a molecule is an enol (vinylic alcohol).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an enol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate through all atoms to find the enol group
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and atom.GetDegree() == 1:
            # Check if the oxygen is connected to a carbon
            carbon_atom = atom.GetNeighbors()[0]
            if carbon_atom.GetSymbol() == 'C' and carbon_atom.GetDegree() == 3:
                # Check if the carbon is part of a double bond (alkene)
                for neighbor in carbon_atom.GetNeighbors():
                    bond = mol.GetBondBetweenAtoms(carbon_atom.GetIdx(), neighbor.GetIdx())
                    if neighbor.GetSymbol() == 'C' and bond and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        # Check if the other carbon in the double bond has at least one hydrogen or is connected to another carbon
                        if any(sub_neighbor.GetSymbol() == 'H' for sub_neighbor in neighbor.GetNeighbors()) or \
                           any(sub_neighbor.GetSymbol() == 'C' for sub_neighbor in neighbor.GetNeighbors()):
                            return True, "Enol group found"
    return False, "No enol group found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33823',
                          'name': 'enol',
                          'definition': 'Alkenols; the term refers '
                                        'specifically to vinylic alcohols, '
                                        "which have the structure HOCR'=CR2. "
                                        'Enols are tautomeric with aldehydes '
                                        "(R' = H) or ketones (R' =/= H).",
                          'parents': ['CHEBI:33822', 'CHEBI:78840']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 19,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 4,
    'precision': 0.95,
    'recall': 0.8260869565217391,
    'f1': 0.8837209302325583,
    'accuracy': None}
"""
Classifies: CHEBI:32878 alkene
"""
"""
Classifies: CHEBI:16336 alkene
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene based on its SMILES string.
    An alkene is an acyclic branched or unbranched hydrocarbon having one carbon-carbon double bond
    and the general formula CnH2n.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly one double bond
    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bonds != 1:
        return False, f"Found {double_bonds} double bonds, must have exactly 1"

    # Check for acyclicity
    if not mol.GetRingInfo().IsAtomacyclic():
        return False, "Molecule is cyclic"

    # Check for only C and H atoms
    atom_types = set(atom.GetAtomicNum() for atom in mol.GetAtoms())
    if atom_types != {1, 6}:
        return False, "Molecule contains atoms other than C and H"

    # Check molecular formula
    formula = rdMolDescriptors.CalcMolFormula(mol)
    n_carbons = formula.count("C")
    n_hydrogens = formula.count("H")
    if n_hydrogens != 2 * n_carbons:
        return False, "Molecular formula does not match CnH2n"

    # Check for branching
    atoms_by_degree = {degree: [atom for atom in mol.GetAtoms() if atom.GetDegree() == degree] for degree in range(1, 5)}
    if len(atoms_by_degree[4]) > 0:
        return True, "Branched alkene"
    elif len(atoms_by_degree[3]) == 2 and len(atoms_by_degree[2]) >= 2:
        return True, "Unbranched alkene"
    else:
        return False, "Invalid alkene structure"

    return True, "Valid alkene molecule"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:16336',
        'name': 'alkene',
        'definition': 'An acyclic branched or unbranched hydrocarbon having one carbon-carbon double bond and the general formula CnH2n. Acyclic branched or unbranched hydrocarbons having more than one double bond are alkadienes, alkatrienes, etc.',
        'parents': ['CHEBI:24703', 'CHEBI:33899']
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 67,
    'num_false_positives': 0,
    'num_true_negatives': 182444,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': 1.0
}
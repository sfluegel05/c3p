"""
Classifies: CHEBI:61384 sulfolipid
"""
"""
Classifies: CHEBI:63505 sulfolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    A sulfolipid is a compound containing a sulfonic acid residue joined by a carbon-sulfur bond to a lipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sulfonic acid group (-SO3H)
    sulfonic_acid_pattern = Chem.MolFromSmarts("[SX4](=[OX1])(=[OX1])([OX2H])")
    if not mol.HasSubstructMatch(sulfonic_acid_pattern):
        return False, "No sulfonic acid group found"

    # Check for a carbon-sulfur bond (C-S) connected to the sulfonic acid group
    carbon_sulfur_bond_pattern = Chem.MolFromSmarts("[CX4]-[SX4](=[OX1])(=[OX1])([OX2H])")
    if not mol.HasSubstructMatch(carbon_sulfur_bond_pattern):
        return False, "No carbon-sulfur bond connected to sulfonic acid group"

    # Check for lipid-like structure (long hydrocarbon chains)
    # We can approximate this by checking for a minimum number of carbons and a certain ratio of carbons to other atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)

    # A lipid typically has a high carbon-to-oxygen ratio
    if c_count < 10 or c_count / (o_count + s_count) < 2:
        return False, "Molecule does not have a lipid-like structure (insufficient carbon content)"

    # Check for long hydrocarbon chains (at least 10 carbons in a row)
    long_chain_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long hydrocarbon chain found"

    # Check molecular weight - sulfolipids typically have a higher molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for sulfolipid"

    return True, "Contains a sulfonic acid group connected to a lipid-like structure via a carbon-sulfur bond"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:63505',
                          'name': 'sulfolipid',
                          'definition': 'A compound containing a sulfonic acid residue joined by a carbon-sulfur bond to a lipid.',
                          'parents': ['CHEBI:63505', 'CHEBI:63505']},
    'config': {   'llm_model_name': 'lbl/claude-sonnet',
                  'f1_threshold': 0.8,
                  'max_attempts': 5,
                  'max_positive_instances': None,
                  'max_positive_to_test': None,
                  'max_negative_to_test': None,
                  'max_positive_in_prompt': 50,
                  'max_negative_in_prompt': 20,
                  'max_instances_in_prompt': 100,
                  'test_proportion': 0.1},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199}
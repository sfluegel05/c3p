"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA (fatty acyl group with chain length greater than C22).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the substructure for CoA
    coa_substructure = Chem.MolFromSmarts('SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12')

    if not mol.HasSubstructMatch(coa_substructure):
        return False, "Molecule does not contain CoA substructure"

    # Remove the CoA substructure to isolate the fatty acyl part
    fatty_acyl_part = Chem.DeleteSubstructs(mol, coa_substructure)

    # Count the number of carbon atoms in the fatty acyl chain
    num_carbons = sum(1 for atom in fatty_acyl_part.GetAtoms() if atom.GetSymbol() == 'C' and atom.GetDegree() < 4)

    if num_carbons > 22:
        return True, f"Fatty acyl chain has {num_carbons} carbon atoms, which is greater than 22"
    else:
        return False, f"Fatty acyl chain has {num_carbons} carbon atoms, which is not greater than 22"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:61910',
                          'name': 'very long-chain fatty acyl-CoA',
                          'definition': 'A fatty acyl-CoA in which the fatty '
                                        'acyl group has a chain length greater '
                                        'than C22.',
                          'parents': ['CHEBI:37554']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 13,
    'num_false_positives': 0,
    'num_true_negatives': 13,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}
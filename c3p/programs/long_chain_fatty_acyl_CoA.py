"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMILES pattern for CoA
    coa_pattern = Chem.MolFromSmarts("C(C(=O)NCCSC(=O))NCCNC(=O)C(O)C(C)(C)COP(O)(=O)OP(O)(=O)OCC1O[C@H](N2C=NC=3C(N)=NC=NC23)C(O)C1OP(O)(O)=O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Molecule does not contain the CoA moiety"

    # Extract the fatty acid part
    fatty_acid_part = None
    for match in mol.GetSubstructMatches(coa_pattern):
        fatty_acid_part = Chem.FragmentOnBonds(mol, [match[0]], addDummies=False)
        break

    if not fatty_acid_part:
        return False, "Could not extract fatty acid part"

    # Check for long-chain fatty acid (C13 to C22)
    carbon_count = 0
    for atom in fatty_acid_part.GetAtoms():
        if atom.GetSymbol() == 'C':
            carbon_count += 1

    if carbon_count < 13 or carbon_count > 22:
        return False, f"Fatty acid chain length is {carbon_count}, which is outside the range C13 to C22"

    return True, "Molecule is a long-chain fatty acyl-CoA"

# Example usage
smiles = "CC/C=C\C/C=C\C/C=C\C/C=C\CCCCC(O)CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1O[C@H](N2C=NC=3C(N)=NC=NC23)C(O)C1OP(=O)(O)O"
result, reason = is_long_chain_fatty_acyl_CoA(smiles)
print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33184',
                          'name': 'long-chain fatty acyl-CoA',
                          'definition': 'A fatty acyl-CoA that results from '
                                        'the formal condensation of the thiol '
                                        'group of coenzyme A with the carboxy '
                                        'group of any long-chain (C13 to C22) '
                                        'fatty acid.',
                          'parents': ['CHEBI:37554']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'False Molecule does not contain the CoA moiety\n',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 20,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}
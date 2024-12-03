"""
Classifies: CHEBI:35748 fatty acid ester
"""
from rdkit import Chem

def is_fatty_acid_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid ester.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    ester_found = False
    fatty_acid_found = False

    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and atom1.GetSymbol() == 'C' and atom2.GetSymbol() == 'O':
            # Check for ester group (C-O-C=O)
            if any(nb.GetSymbol() == 'C' and any(nnb.GetSymbol() == 'O' and nnb.GetIdx() != atom2.GetIdx() for nnb in nb.GetNeighbors()) for nb in atom1.GetNeighbors()):
                ester_found = True
                break

    if ester_found:
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'C' and atom.GetDegree() > 1:
                # Check if it is part of a fatty acid chain (long carbon chain)
                chain_length = 0
                visited = set()
                stack = [(atom, 0)]
                while stack:
                    current_atom, length = stack.pop()
                    if current_atom.GetIdx() not in visited:
                        visited.add(current_atom.GetIdx())
                        if current_atom.GetSymbol() == 'C':
                            chain_length = max(chain_length, length + 1)
                        for neighbor in current_atom.GetNeighbors():
                            if neighbor.GetSymbol() == 'C' and neighbor.GetIdx() not in visited:
                                stack.append((neighbor, length + 1))
                if chain_length >= 8:  # Assuming a minimum chain length for fatty acids
                    fatty_acid_found = True
                    break

    if ester_found and fatty_acid_found:
        return True, "Molecule is a fatty acid ester"
    elif ester_found:
        return False, "Ester group found but no long carbon chain"
    else:
        return False, "No ester group found"

# Example usage
smiles = "OCC(COC(CCCCCCCCCCCCC/C=C\\CCCCCCCC)=O)O"
print(is_fatty_acid_ester(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35748',
                          'name': 'fatty acid ester',
                          'definition': 'A carboxylic ester in which the '
                                        'carboxylic acid component can be any '
                                        'fatty acid.',
                          'parents': ['CHEBI:33308', 'CHEBI:61697']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(True, 'Molecule is a fatty acid ester')\n",
    'num_true_positives': 48,
    'num_false_positives': 4,
    'num_true_negatives': 16,
    'num_false_negatives': 85,
    'precision': 0.9230769230769231,
    'recall': 0.3609022556390977,
    'f1': 0.518918918918919,
    'accuracy': None}
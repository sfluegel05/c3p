"""
Classifies: CHEBI:16337 phosphatidic acid
"""
"""
Classifies: phosphatidic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    A phosphatidic acid is a derivative of glycerol in which one hydroxy group is esterified with phosphoric acid
    and the other two are esterified with fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphate group
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)[O]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Look for glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("[CH2][CH][CH2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for ester groups connected to glycerol carbons
    ester_pattern = Chem.MolFromSmarts("C(=O)O[CH]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups connected to glycerol carbons, need at least 2"

    # Check that esters are connected to fatty acid chains (long carbon chains)
    long_chain_count = 0
    for match in ester_matches:
        carbonyl_c_idx = match[0]
        # Get the carbonyl carbon atom
        c_atom = mol.GetAtomWithIdx(carbonyl_c_idx)
        # Get the carbon chain length starting from the carbonyl carbon
        chain_length = count_chain_length(c_atom, mol)
        if chain_length >= 8:
            long_chain_count += 1

    if long_chain_count < 2:
        return False, f"Found {long_chain_count} fatty acid chains of sufficient length, need at least 2"

    return True, "Contains glycerol backbone with phosphate group and two fatty acid chains attached via ester bonds"

# Helper function to count chain length
def count_chain_length(c_atom, mol):
    visited = set()
    queue = [(c_atom, 0)]
    max_chain_length = 0
    while queue:
        atom, length = queue.pop(0)
        atom_idx = atom.GetIdx()
        if atom_idx in visited:
            continue
        visited.add(atom_idx)
        if atom.GetAtomicNum() == 6:
            length += 1
            if length > max_chain_length:
                max_chain_length = length
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                    queue.append((neighbor, length))
    return max_chain_length

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:49105',
                              'name': 'phosphatidic acid',
                              'definition': 'A derivative of glycerol in which one hydroxy group, commonly but not necessarily primary, is esterified with phosphoric acid and the other two are esterified with fatty acids.',
                              'parents': []},
        'config': {   'llm_model_name': 'lbl/claude-sonnet',
                      'f1_threshold': 0.8,
                      'max_attempts': 5},
        'message': None,
        'attempt': 0,
        'success': True,
        'best': True,
        'error': '',
        'stdout': None}
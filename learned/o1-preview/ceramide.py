"""
Classifies: CHEBI:17761 ceramide
"""
"""
Classifies: CHEBI:17761 ceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def get_chain_length_from_atom(mol, start_atom, exclude_atoms=[]):
    visited = set(exclude_atoms)
    max_length = [0]

    def dfs(atom, length):
        atom_idx = atom.GetIdx()
        if atom_idx in visited:
            return
        visited.add(atom_idx)
        if atom.GetSymbol() != 'C':
            return
        length += 1
        if length > max_length[0]:
            max_length[0] = length
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in visited:
                dfs(neighbor, length)

    dfs(start_atom, 0)
    return max_length[0]

def get_sphingoid_base_info(mol, start_atom, exclude_atoms=[]):
    visited = set(exclude_atoms)
    num_hydroxy = 0
    num_amine = 0
    chain_length = 0

    def dfs(atom):
        nonlocal num_hydroxy, num_amine, chain_length
        atom_idx = atom.GetIdx()
        if atom_idx in visited:
            return
        visited.add(atom_idx)
        symbol = atom.GetSymbol()
        if symbol == 'C':
            chain_length += 1
        elif symbol == 'O':
            num_hydroxy += 1
        elif symbol == 'N' and atom_idx != start_atom.GetIdx():
            num_amine += 1
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in visited:
                dfs(neighbor)

    dfs(start_atom)
    return chain_length, num_hydroxy, num_amine

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.
    Ceramides are sphingoid bases acylated via an amide linkage to a fatty acid chain.
    The sphingoid base typically has two hydroxyl groups and an amino group, forming a long-chain amino alcohol.
    The fatty acid chain is typically saturated or monounsaturated, with chain lengths from 14 to 26 carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a ceramide, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find amide bonds
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)

    if not amide_matches:
        return False, "No amide bond found"

    # Assume that ceramide has only one amide bond
    if len(amide_matches) != 1:
        return False, f"Expected one amide bond, found {len(amide_matches)}"

    # Get the amide bond atoms
    amide_match = amide_matches[0]
    carbonyl_c_idx = amide_match[0]
    nitrogen_idx = amide_match[2]

    # Get the fatty acid chain length (connected to the carbonyl carbon)
    carbonyl_c_atom = mol.GetAtomWithIdx(carbonyl_c_idx)
    exclude_atoms = [amide_match[1], nitrogen_idx]
    fatty_acid_chain_length = get_chain_length_from_atom(mol, carbonyl_c_atom, exclude_atoms=exclude_atoms)

    # Check if fatty acid chain length is between 14 and 26 carbons
    if not (14 <= fatty_acid_chain_length <= 26):
        return False, f"Fatty acid chain length is {fatty_acid_chain_length}, expected between 14 and 26 carbons"

    # Get the sphingoid base info (connected to the nitrogen)
    nitrogen_atom = mol.GetAtomWithIdx(nitrogen_idx)
    exclude_atoms = [carbonyl_c_idx]
    sphingoid_chain_length, num_hydroxy, num_amine = get_sphingoid_base_info(mol, nitrogen_atom, exclude_atoms=exclude_atoms)

    if num_hydroxy < 2:
        return False, f"Sphingoid base has {num_hydroxy} hydroxyl groups, expected at least 2"

    if sphingoid_chain_length < 12:
        return False, f"Sphingoid base chain length is {sphingoid_chain_length}, expected at least 12 carbons"

    return True, "Molecule matches ceramide definition: sphingoid base with amide-linked fatty acid chain"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:17761',
        'name': 'ceramide',
        'definition': "Ceramides (N-acyl-sphingoid bases) are a major subclass of sphingoid base derivatives with an amide-linked fatty acid. The fatty acids are typically saturated or monounsaturated with chain lengths from 14 to 26 carbon atoms; the presence of a hydroxyl group on carbon 2 is fairly common. Ceramides are generally precursors of more complex sphingolipids.",
        'parents': ['CHEBI:33536']
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
    # Statistical measures such as 'num_true_positives', 'precision', 'recall', etc., would be filled after testing the classifier
}
"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
"""
Classifies: CHEBI:52634 N-acylsphinganine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.
    An N-acylsphinganine is a ceramide consisting of sphinganine in which one of the amino hydrogens is substituted by a fatty acyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphinganine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define sphinganine backbone pattern
    # Sphinganine backbone: C1-C2(NH)-C3(OH)-[long chain]-C(OH)
    sphinganine_pattern = Chem.MolFromSmarts("""
        [C;H2][C@H](N)[C@@H](O)[C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C](O)
    """)

    if not mol.HasSubstructMatch(sphinganine_pattern):
        return False, "No sphinganine backbone found"

    # Define N-acyl group pattern (amide linkage with fatty acyl chain)
    # Amide bond: [N]-C(=O)-[C;H2][C;H2][C;H2]... (long chain)
    acylamide_pattern = Chem.MolFromSmarts("""
        [NX3][C](=O)[C;H2][C;H2][C;H2][C;H2][C;H2]
    """)

    if not mol.HasSubstructMatch(acylamide_pattern):
        return False, "No N-acyl group (amide bond with fatty acyl chain) found"

    # Ensure that the nitrogen in the sphinganine backbone is the one acylated
    sphinganine_matches = mol.GetSubstructMatches(sphinganine_pattern)
    acylamide_matches = mol.GetSubstructMatches(acylamide_pattern)

    sphinganine_nitrogens = set()
    for match in sphinganine_matches:
        sphinganine_nitrogen_idx = match[2]  # The nitrogen atom in sphinganine pattern
        sphinganine_nitrogens.add(sphinganine_nitrogen_idx)

    acylated_nitrogens = set()
    for match in acylamide_matches:
        acyl_nitrogen_idx = match[0]  # The nitrogen atom in acylamide pattern
        acylated_nitrogens.add(acyl_nitrogen_idx)

    # Check if the sphinganine nitrogen is acylated
    common_nitrogens = sphinganine_nitrogens & acylated_nitrogens
    if not common_nitrogens:
        return False, "The sphinganine nitrogen is not acylated"

    # Optional: Check that the acyl chain is sufficiently long (e.g., at least 12 carbons)
    acyl_chain_lengths = []
    for match in acylamide_matches:
        carbonyl_carbon = mol.GetAtomWithIdx(match[1])
        acyl_chain_atom = mol.GetAtomWithIdx(match[2])

        # Traverse the acyl chain
        acyl_chain_length = 0
        visited = set()
        stack = [acyl_chain_atom]
        while stack:
            atom = stack.pop()
            if atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())
            if atom.GetAtomicNum() == 6:  # Carbon
                acyl_chain_length += 1
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                        stack.append(neighbor)
        acyl_chain_lengths.append(acyl_chain_length)

    if not any(length >= 12 for length in acyl_chain_lengths):
        return False, "Acyl chain is too short to be a fatty acyl group"

    return True, "Contains sphinganine backbone with N-acyl group forming an amide bond"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:52634',
        'name': 'N-acylsphinganine',
        'definition': "A ceramide consisting of sphinganine in which one of the amino hydrogens is substituted by a fatty acyl group.",
        'parents': ['CHEBI:33598', 'CHEBI:33709']
    },
    'config': {}
}
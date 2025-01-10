"""
Classifies: CHEBI:64583 sphingomyelin
"""
"""
Classifies: CHEBI:15837 sphingomyelin
"""

from rdkit import Chem

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    A sphingomyelin is characterized by:
    - A sphingoid base backbone: a long-chain amino alcohol (variable length, usually 18 carbons) with two hydroxyl groups and one amino group.
    - An amide linkage between the amino group of the sphingoid base and a fatty acid.
    - A phosphocholine group attached to the terminal hydroxyl group via a phosphodiester linkage.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a sphingomyelin, False otherwise.
        str: Reason for classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure molecule is not empty
    if mol.GetNumAtoms() == 0:
        return False, "Empty molecule"

    # Define SMARTS patterns

    # 1. Sphingoid base backbone pattern
    # Long hydrocarbon chain with two hydroxyl groups and one amino group
    # Allow for variable chain length, unsaturation, and possible methyl branches
    sphingoid_pattern = Chem.MolFromSmarts("""
        [$([C;H2,H1](O))]                # Primary or secondary carbon with hydroxyl group (start of sphingoid base)
        [C;!R;D2,D3][$([C;!R])]*
        [C;!R;D3](N)                     # Carbon connected to amino group
        [$([C;!R;D2,D3][$([C;!R])]*)]*
        [$([C;H2,H1](O))]                # Primary or secondary carbon with hydroxyl group (end of sphingoid base)
    """.replace('\n', '').replace(' ', ''))

    if not mol.HasSubstructMatch(sphingoid_pattern):
        return False, "No sphingoid base backbone found"

    # Get the atoms involved in the sphingoid base backbone
    sphingoid_match = mol.GetSubstructMatch(sphingoid_pattern)
    sphingoid_atoms = set(sphingoid_match)

    # 2. Amide linkage pattern (amide bond connected to sphingoid amino group)
    amide_pattern = Chem.MolFromSmarts("""
        [N;!R]C(=O)[C]                    # Amide nitrogen connected to carbonyl carbon
    """.replace('\n', '').replace(' ', ''))

    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide linkage with fatty acid found"

    # Check if amide nitrogen is part of sphingoid base
    amide_nitrogen_idx = None
    for match in amide_matches:
        nitrogen_idx = match[0]
        if nitrogen_idx in sphingoid_atoms:
            amide_nitrogen_idx = nitrogen_idx
            break
    if amide_nitrogen_idx is None:
        return False, "Amide linkage not connected to sphingoid amino group"

    # 3. Phosphocholine group pattern
    phosphocholine_pattern = Chem.MolFromSmarts("""
        [O]-P(=O)([O-])OCC[N+](C)(C)C     # Phosphocholine group
    """.replace('\n', '').replace(' ', ''))

    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group found"

    # Check if phosphocholine group is connected to sphingoid terminal hydroxyl
    phosphocholine_matches = mol.GetSubstructMatches(phosphocholine_pattern)
    connected = False
    for match in phosphocholine_matches:
        oxygen_idx = match[0]  # The oxygen connected to the phosphorus
        oxygen_atom = mol.GetAtomWithIdx(oxygen_idx)
        neighbors = oxygen_atom.GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetIdx() in sphingoid_atoms and neighbor.GetAtomicNum() == 6:
                connected = True
                break
        if connected:
            break
    if not connected:
        return False, "Phosphocholine group not connected to sphingoid terminal hydroxyl"

    return True, "Molecule is a sphingomyelin with correct structural features"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:15837',
        'name': 'sphingomyelin',
        'definition': 'Any of a class of phospholipids in which the amino group of a sphingoid base is in amide linkage with one of several fatty acids, while the terminal hydroxy group of the sphingoid base is esterified to phosphorylcholine.',
        'parents': ['CHEBI:64713', 'CHEBI:65409']
    }
}
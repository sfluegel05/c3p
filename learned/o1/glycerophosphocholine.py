"""
Classifies: CHEBI:36313 glycerophosphocholine
"""
"""
Classifies: glycerophosphocholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycerophosphocholine(smiles: str):
    """
    Determines if a molecule is a glycerophosphocholine based on its SMILES string.
    A glycerophosphocholine has a glycerol backbone with fatty acid chains attached
    at positions 1 and 2 via ester or ether bonds, and a phosphocholine group at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphocholine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns with atom mapping
    glycerol_pattern = Chem.MolFromSmarts("""
    [C:1]-[C@H:2]-[C:3]
    """)
    phosphocholine_pattern = Chem.MolFromSmarts("""
    [O:4]-[P:5](=O)([O-:6])-[O:7]-[C:8][C:9][N+:10]([C:11])([C:12])[C:13]
    """)
    ester_or_ether_pattern = Chem.MolFromSmarts("""
    [C:14]-[O:15]-[C:16]
    """)

    # Find matches for glycerol backbone
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "No glycerol backbone found"

    # For each glycerol match, attempt to find attached groups
    for glycerol_match in glycerol_matches:
        c1_idx = glycerol_match[0]
        c2_idx = glycerol_match[1]
        c3_idx = glycerol_match[2]

        # Check for ester or ether groups at positions 1 and 2
        fatty_acid_chains = 0
        for c_idx in [c1_idx, c2_idx]:
            found_chain = False
            for bond in mol.GetAtomWithIdx(c_idx).GetBonds():
                neighbor_atom = bond.GetOtherAtom(mol.GetAtomWithIdx(c_idx))
                neighbor_idx = neighbor_atom.GetIdx()
                if neighbor_atom.GetSymbol() == 'O':
                    # Check if this oxygen is part of an ester or ether
                    for pattern in [Chem.MolFromSmarts("[O]-C=O"), Chem.MolFromSmarts("[O]-[C]")]:
                        submol = Chem.PathToSubmol(mol, [c_idx, neighbor_idx])
                        if submol.HasSubstructMatch(pattern):
                            # Check for long carbon chain attached
                            attached_atom = None
                            if neighbor_atom.GetDegree() > 1:
                                for neighbor in neighbor_atom.GetNeighbors():
                                    if neighbor.GetIdx() != c_idx:
                                        attached_atom = neighbor
                                        break
                            if attached_atom and attached_atom.GetSymbol() == 'C':
                                chain_length = CountCarbonChain(mol, attached_atom.GetIdx(), exclude_atoms=[c_idx, neighbor_idx])
                                if chain_length >= 4:
                                    fatty_acid_chains += 1
                                    found_chain = True
                                    break
                if found_chain:
                    break

        if fatty_acid_chains < 2:
            continue  # Not enough fatty acid chains attached at positions 1 and 2

        # Check for phosphocholine group at position 3
        phosphate_found = False
        for bond in mol.GetAtomWithIdx(c3_idx).GetBonds():
            neighbor_atom = bond.GetOtherAtom(mol.GetAtomWithIdx(c3_idx))
            if neighbor_atom.GetSymbol() == 'O':
                # Check if this oxygen is connected to a phosphate group
                phosphate_atom = None
                for neighbor in neighbor_atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'P':
                        phosphate_atom = neighbor
                        break
                if phosphate_atom:
                    # Check if phosphate group matches phosphocholine pattern
                    match = mol.HasSubstructMatch(phosphocholine_pattern)
                    if match:
                        phosphate_found = True
                        break
        if not phosphate_found:
            continue  # No phosphocholine group attached at position 3

        return True, "Molecule is a glycerophosphocholine with glycerol backbone, fatty acid chains, and phosphocholine group"

    return False, "Does not match glycerophosphocholine structure"

def CountCarbonChain(mol, start_idx, exclude_atoms=None):
    """
    Count the number of carbons in a chain starting from a given atom index.

    Args:
        mol (Chem.Mol): RDKit molecule object
        start_idx (int): Atom index to start counting from
        exclude_atoms (list): List of atom indices to exclude from counting

    Returns:
        int: Number of carbons in the chain
    """
    if exclude_atoms is None:
        exclude_atoms = []
    visited = set(exclude_atoms)
    count = 0
    stack = [start_idx]
    while stack:
        idx = stack.pop()
        if idx in visited:
            continue
        visited.add(idx)
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetSymbol() == 'C':
            count += 1
            for neighbor in atom.GetNeighbors():
                stack.append(neighbor.GetIdx())
    return count

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'glycerophosphocholine',
        'definition': 'The glycerol phosphate ester of a phosphocholine. A nutrient with many different roles in human health.',
        'parents': []
    }
}
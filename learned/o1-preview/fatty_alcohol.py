"""
Classifies: CHEBI:24026 fatty alcohol
"""
"""
Classifies: fatty alcohol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.
    A fatty alcohol is an aliphatic alcohol consisting of a chain of 3 to greater than 27 carbon atoms.
    Fatty alcohols may be saturated or unsaturated and may be branched or unbranched.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty alcohol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of hydroxyl group (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group (-OH) found"

    # Check that the molecule is aliphatic (no aromatic rings)
    if mol.GetRingInfo().NumAromaticRings() > 0:
        return False, "Contains aromatic rings, not aliphatic"

    # Get all carbon atoms
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]

    if len(carbon_atoms) < 3:
        return False, f"Contains only {len(carbon_atoms)} carbon atoms, minimum is 3"

    # Find the longest aliphatic carbon chain
    # First, sanitize molecule to kekulize any aromatic features
    Chem.Kekulize(mol, clearAromaticFlags=True)
    paths = Chem.FindAllPathsOfLengthN(mol, 1, False)  # Initialize paths
    max_chain_length = 0
    for bond in mol.GetBonds():
        # Consider bonds between carbon atoms
        if bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6:
            visited_atoms = set()
            stack = [(bond.GetBeginAtomIdx(), 1)]
            while stack:
                current_atom_idx, length = stack.pop()
                visited_atoms.add(current_atom_idx)
                max_chain_length = max(max_chain_length, length)
                if length > 100:  # To prevent infinite loops in cyclic structures
                    break
                atom = mol.GetAtomWithIdx(current_atom_idx)
                for neighbor in atom.GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor.GetAtomicNum() == 6 and neighbor_idx not in visited_atoms:
                        stack.append((neighbor_idx, length + 1))

    if max_chain_length < 3:
        return False, f"Longest carbon chain is {max_chain_length}, minimum is 3"

    # All conditions met
    return True, "Molecule is an aliphatic alcohol with at least 3 carbons"

__metadata__ = {
    'chemical_class': {
        'name': 'fatty alcohol',
        'definition': 'An aliphatic alcohol consisting of a chain of 3 to greater than 27 carbon atoms. Fatty alcohols may be saturated or unsaturated and may be branched or unbranched.',
    }
}
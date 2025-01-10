"""
Classifies: CHEBI:35785 sphingoid
"""
"""
Classifies: sphingoid
"""
from rdkit import Chem

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    A sphingoid is defined as sphinganine, its homologs and stereoisomers,
    and the hydroxy and unsaturated derivatives of these compounds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find the longest aliphatic carbon chain (sphingoid backbone)
    def get_longest_chain(mol):
        chains = []
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6 and not atom.IsInRing():
                visited = set()
                stack = [(atom, [atom])]
                while stack:
                    current_atom, path = stack.pop()
                    visited.add(current_atom.GetIdx())
                    neighbors = [nbr for nbr in current_atom.GetNeighbors()
                                 if nbr.GetAtomicNum() == 6 and not nbr.IsInRing() and nbr.GetIdx() not in visited]
                    if neighbors:
                        for nbr in neighbors:
                            stack.append((nbr, path + [nbr]))
                    else:
                        chains.append(path)
        if not chains:
            return []
        # Return the longest chain (list of atoms)
        longest_chain = max(chains, key=lambda x: len(x))
        return longest_chain

    backbone = get_longest_chain(mol)
    if len(backbone) < 12:
        return False, f"Aliphatic chain too short ({len(backbone)} carbons)"

    # Map positions: C1, C2, C3,...
    # C1 is the first carbon in the chain
    c1 = backbone[0]
    c2 = backbone[1] if len(backbone) > 1 else None
    c3 = backbone[2] if len(backbone) > 2 else None

    # Check for hydroxyl group on C1
    has_oh_c1 = False
    for neighbor in c1.GetNeighbors():
        if neighbor.GetAtomicNum() == 8:
            has_oh_c1 = True
            break

    # Check for amino group on C2 (including protonated forms)
    has_n_c2 = False
    if c2:
        for neighbor in c2.GetNeighbors():
            if neighbor.GetAtomicNum() == 7:
                has_n_c2 = True
                break

    # Check for hydroxyl group on C3
    has_oh_c3 = False
    if c3:
        for neighbor in c3.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:
                has_oh_c3 = True
                break

    if not has_n_c2:
        return False, "No amino group attached to C2 of the backbone"

    if not (has_oh_c1 or has_oh_c3):
        return False, "No hydroxyl groups attached to C1 or C3 of the backbone"

    # Allow for unsaturation and additional hydroxyl groups (derivatives)
    # Check that the molecule is not a peptide or contains other functional groups that exclude it from being a sphingoid
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    if mol.HasSubstructMatch(peptide_bond_pattern):
        return False, "Contains peptide bonds, not a sphingoid"

    # If all checks pass, classify as sphingoid
    return True, "Molecule matches sphingoid structural features"

__metadata__ = {
    'chemical_class': {
        'name': 'sphingoid',
        'definition': 'Sphinganine, its homologs and stereoisomers, and the hydroxy and unsaturated derivatives of these compounds.'
    }
}
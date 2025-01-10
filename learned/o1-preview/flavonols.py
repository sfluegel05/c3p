"""
Classifies: CHEBI:28802 flavonols
"""
"""
Classifies: flavonols
"""

from rdkit import Chem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    A flavonol is defined as any hydroxyflavone in which the hydrogen at position 3
    of the heterocyclic ring is replaced by a hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for flavone core (2-phenylchromen-4-one)
    flavone_smarts = 'O=C1C=CC2=CC=CC=C2O1'
    flavone_pattern = Chem.MolFromSmarts(flavone_smarts)
    if flavone_pattern is None:
        return False, "Invalid flavone SMARTS pattern"

    # Find matches of flavone core
    matches = mol.GetSubstructMatches(flavone_pattern)
    if not matches:
        return False, "Flavone core not found"

    # For each match of the flavone core
    for match in matches:
        # Build atom mapping from query to target atoms
        atom_map = {flavone_atom_idx: mol_atom_idx for flavone_atom_idx, mol_atom_idx in enumerate(match)}

        # Identify the atom at position 3 in the flavone core
        # In the flavone SMARTS, atom indices are as follows:
        # Index 0: O (carbonyl oxygen at position 4)
        # Index 1: C (carbon at position 4)
        # Index 2: C (carbon at position 3)
        # Index 3: C (carbon at position 2)
        # Index 4-9: Phenyl ring attached at position 2

        # Get the atom index for position 3
        position_3_atom_idx = atom_map[2]  # Atom at position 3 in the flavone core
        position_3_atom = mol.GetAtomWithIdx(position_3_atom_idx)

        # Check if the atom at position 3 has a hydroxy group attached
        has_3_hydroxy = False
        for neighbor in position_3_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen atom
                # Check if oxygen is connected via single bond
                bond = mol.GetBondBetweenAtoms(position_3_atom_idx, neighbor.GetIdx())
                if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    # Check if oxygen is part of a hydroxy group (has hydrogen attached)
                    num_h = neighbor.GetTotalNumHs()
                    if num_h == 1:
                        has_3_hydroxy = True
                        break
        if has_3_hydroxy:
            return True, "Contains flavone core with hydroxy group at position 3"

    return False, "Flavone core found but no hydroxy group at position 3"
"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
"""
Classifies: CHEBI:24996 octadecadienoic acid
"""

from rdkit import Chem

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.
    An octadecadienoic acid is any straight-chain C18 polyunsaturated fatty acid
    having two C=C double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octadecadienoic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group and get the carboxyl carbon atom
    carboxylic_acid = Chem.MolFromSmarts("C(=O)[O;H1,H0-]")
    matches = mol.GetSubstructMatch(carboxylic_acid)
    if not matches:
        return False, "No carboxylic acid group found"
    carboxyl_c_idx = matches[0]  # Index of the carboxyl carbon atom
    carboxyl_c = mol.GetAtomWithIdx(carboxyl_c_idx)

    # Ensure the molecule is acyclic (no rings)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings; not a straight-chain fatty acid"

    # Start traversal from the carbon next to the carboxyl carbon
    neighbors = [nbr for nbr in carboxyl_c.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(neighbors) != 1:
        return False, "Carboxyl carbon should be connected to exactly one carbon atom"
    current_atom = neighbors[0]
    previous_atom = carboxyl_c

    chain_length = 2  # Starts from carboxyl carbon and next carbon
    double_bonds = 0
    is_branched = False

    # Check bond order between carboxyl carbon and first carbon
    bond = mol.GetBondBetweenAtoms(carboxyl_c_idx, current_atom.GetIdx())
    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
        double_bonds += 1

    visited_atoms = set()
    visited_atoms.add(carboxyl_c_idx)
    visited_atoms.add(current_atom.GetIdx())

    while True:
        # Get carbon neighbors excluding the previous atom
        neighbors = [nbr for nbr in current_atom.GetNeighbors()
                     if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != previous_atom.GetIdx()]
        if len(neighbors) == 0:
            # Reached the end of the chain
            break
        if len(neighbors) > 1:
            # Branching detected
            is_branched = True
            break  # Not a straight-chain fatty acid

        next_atom = neighbors[0]

        # Check for branching at the next atom
        next_neighbors = [nbr for nbr in next_atom.GetNeighbors()
                          if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != current_atom.GetIdx()]
        if len(next_neighbors) > 1:
            # Branching detected at the next atom
            is_branched = True
            break

        # Check bond order between current atom and next atom
        bond = mol.GetBondBetweenAtoms(current_atom.GetIdx(), next_atom.GetIdx())
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            double_bonds += 1

        # Update traversal
        previous_atom = current_atom
        current_atom = next_atom
        visited_atoms.add(current_atom.GetIdx())
        chain_length += 1

    if is_branched:
        return False, "Molecule is branched; not a straight-chain fatty acid"

    if chain_length != 18:
        return False, f"Chain length is {chain_length}, expected 18 carbons"

    if double_bonds != 2:
        return False, f"Number of C=C double bonds along the chain is {double_bonds}, expected 2"

    return True, "Molecule has a straight-chain C18 fatty acid backbone with two C=C double bonds"
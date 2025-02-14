"""
Classifies: CHEBI:26244 prenols
"""
"""
Classifies: prenols
"""
from rdkit import Chem

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    A prenol is any alcohol possessing the general formula H-[CH2C(Me)=CHCH2]nOH,
    in which the carbon skeleton is composed of one or more isoprene units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for primary alcohol group (-CH2OH)
    OH_pattern = Chem.MolFromSmarts('[C;H2][OX2H]')
    OH_matches = mol.GetSubstructMatches(OH_pattern)
    if not OH_matches:
        return False, "No primary alcohol group (-CH2OH) found"

    # Assume primary alcohol; get carbon attached to OH
    OH_atom_idx = OH_matches[0][1]
    OH_atom = mol.GetAtomWithIdx(OH_atom_idx)
    carbon_atom = mol.GetAtomWithIdx(OH_matches[0][0])

    visited = set()
    chain_atoms = []

    # Function to recursively traverse the main chain
    def traverse_chain(atom, prev_atom_idx, position):
        visited.add(atom.GetIdx())
        chain_atoms.append((atom, position))
        neighbors = [a for a in atom.GetNeighbors() if a.GetIdx() != prev_atom_idx]

        # Exclude hydrogens
        neighbors = [a for a in neighbors if a.GetAtomicNum() != 1]

        main_chain_neighbor = None
        side_chains = []
        for neighbor in neighbors:
            if neighbor.GetIdx() in visited:
                continue
            if neighbor.GetAtomicNum() != 6:
                return False  # Non-carbon atom in chain
            if neighbor.GetDegree() == 1:
                # Methyl group (CH3), acceptable as side chain
                side_chains.append(neighbor)
                continue
            elif main_chain_neighbor is None:
                main_chain_neighbor = neighbor
            else:
                return False  # More than one main chain neighbor, invalid branch

        # Check for side chains at correct positions
        if position % 5 == 2:
            if len(side_chains) != 1:
                return False  # Missing methyl group at position 2 of isoprene unit
        else:
            if len(side_chains) != 0:
                return False  # Unexpected branch at this position

        if main_chain_neighbor is None:
            # End of chain
            return True

        # Recursively traverse the next atom in the main chain
        return traverse_chain(main_chain_neighbor, atom.GetIdx(), position + 1)

    # Start traversing from the carbon attached to OH (position 0)
    is_linear = traverse_chain(carbon_atom, OH_atom_idx, position=0)
    if not is_linear:
        return False, "Chain is branched incorrectly or contains non-carbon atoms"

    # Ensure the chain length corresponds to complete isoprene units
    chain_length = len(chain_atoms)
    n_units = chain_length // 5
    if chain_length % 5 != 0 or n_units < 1:
        return False, f"Chain length ({chain_length}) does not correspond to complete isoprene units"

    # Check double bonds and methyl groups in each isoprene unit
    for i in range(n_units):
        unit_atoms = [atom for atom, pos in chain_atoms[i*5:(i+1)*5]]

        # Positions in the isoprene unit
        pos0 = unit_atoms[0]  # CH2
        pos1 = unit_atoms[1]  # C with double bond to next
        pos2 = unit_atoms[2]  # C with methyl branch
        pos3 = unit_atoms[3]  # C with double bond to next
        pos4 = unit_atoms[4]  # CH2

        # Check for double bond between pos0 and pos1
        bond01 = mol.GetBondBetweenAtoms(pos0.GetIdx(), pos1.GetIdx())
        if bond01.GetBondType() != Chem.rdchem.BondType.SINGLE:
            return False, f"Expected single bond between positions {i*5} and {i*5 + 1}"

        # Check for double bond between pos1 and pos2
        bond12 = mol.GetBondBetweenAtoms(pos1.GetIdx(), pos2.GetIdx())
        if bond12.GetBondType() != Chem.rdchem.BondType.DOUBLE:
            return False, f"Expected double bond between positions {i*5 + 1} and {i*5 + 2}"

        # Check for single bond between pos2 and pos3
        bond23 = mol.GetBondBetweenAtoms(pos2.GetIdx(), pos3.GetIdx())
        if bond23.GetBondType() != Chem.rdchem.BondType.SINGLE:
            return False, f"Expected single bond between positions {i*5 + 2} and {i*5 + 3}"

        # Check for double bond between pos3 and pos4
        bond34 = mol.GetBondBetweenAtoms(pos3.GetIdx(), pos4.GetIdx())
        if bond34.GetBondType() != Chem.rdchem.BondType.DOUBLE:
            return False, f"Expected double bond between positions {i*5 + 3} and {i*5 + 4}"

    return True, "Molecule is a prenol with linear chain of isoprene units ending with -OH group"
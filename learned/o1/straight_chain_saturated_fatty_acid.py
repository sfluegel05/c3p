"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid based on its SMILES string.
    A straight-chain saturated fatty acid is a carboxylic acid with a linear saturated alkyl chain and no side chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a straight-chain saturated fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly one carboxylic acid group (-C(=O)O[H])
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O;H1]')
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 1:
        return False, f"Found {len(carboxylic_acid_matches)} carboxylic acid groups, need exactly 1"

    # Get the carboxyl carbon atom index
    carboxyl_carbon_idx = carboxylic_acid_matches[0][0]

    # Check that the rest of the molecule is a linear, saturated hydrocarbon chain
    # Exclude the carboxyl group oxygens from the traversal
    atom_queue = [carboxyl_carbon_idx]
    visited_atoms = set()
    branching = False
    unsaturated = False

    while atom_queue:
        current_atom_idx = atom_queue.pop()
        if current_atom_idx in visited_atoms:
            continue
        visited_atoms.add(current_atom_idx)
        atom = mol.GetAtomWithIdx(current_atom_idx)

        if atom.GetAtomicNum() == 6:  # Carbon atom
            neighbor_carbons = []
            for bond in atom.GetBonds():
                neighbor_atom = bond.GetOtherAtom(atom)
                neighbor_idx = neighbor_atom.GetIdx()
                if neighbor_atom.GetAtomicNum() == 6:
                    neighbor_carbons.append(neighbor_idx)
                    if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                        unsaturated = True
                elif neighbor_atom.GetAtomicNum() == 8:
                    # Skip carboxyl oxygens
                    continue
                elif neighbor_atom.GetAtomicNum() != 1:
                    # Found a heteroatom (other than oxygen in carboxyl group), molecule is branched
                    branching = True

            # Check for branching (more than 2 neighbor carbons except for terminal carbons)
            if len(neighbor_carbons) > 2:
                branching = True
            elif len(neighbor_carbons) == 2 or (current_atom_idx == carboxyl_carbon_idx and len(neighbor_carbons) == 1):
                # Not a terminal carbon, continue traversal
                atom_queue.extend(neighbor_carbons)
            elif len(neighbor_carbons) == 1:
                # Possible terminal carbon
                atom_queue.extend(neighbor_carbons)
            else:
                # Unexpected number of carbon neighbors
                branching = True

        else:
            # Non-carbon atom (excluding carboxyl oxygens), molecule is branched or contains heteroatoms
            if atom.GetAtomicNum() != 8 or not atom.IsInRing():
                branching = True

    if branching:
        return False, "Molecule is branched or contains heteroatoms"

    if unsaturated:
        return False, "Molecule is unsaturated (contains double or triple bonds)"

    # Check that all atoms have been visited (no disconnected fragments)
    if len(visited_atoms) != mol.GetNumAtoms():
        return False, "Molecule contains disconnected fragments"

    # Check that chain length is at least 4 carbons
    num_carbons = sum(1 for idx in visited_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    if num_carbons < 4:
        return False, f"Chain length {num_carbons} is too short for a fatty acid"

    return True, "Molecule is a straight-chain saturated fatty acid"
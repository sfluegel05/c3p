"""
Classifies: CHEBI:73702 wax
"""
"""
Classifies: wax
"""
from rdkit import Chem

def is_wax(smiles: str):
    """
    Determines if a molecule is a wax based on its SMILES string.
    Waxes are esters formed from long-chain fatty acids and long-chain alcohols, typically with
    carbon chains of 14 carbons or more.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find ester functional groups
    ester_pattern = Chem.MolFromSmarts("[C:1](=O)[O:2][C]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester functional group found"

    # For each ester group, check the length of the carbon chains on both sides
    for match in ester_matches:
        c_atom_idx = match[0]  # Carbonyl carbon atom index
        o_atom_idx = match[1]  # Single-bonded oxygen atom index

        # Break the bond between carbonyl carbon and single-bonded oxygen
        bond = mol.GetBondBetweenAtoms(c_atom_idx, o_atom_idx)
        if bond is None:
            continue  # Bond not found, skip to next match
        bond_idx = bond.GetIdx()

        # Break the bond to create fragments
        mol_frag = Chem.FragmentOnBonds(mol, [bond_idx], addDummies=True)
        fragments = Chem.GetMolFrags(mol_frag, asMols=True, sanitizeFrags=True)

        if len(fragments) != 2:
            continue  # Should result in exactly two fragments, else skip

        # Initialize chain lengths
        acyl_chain_length = 0
        alkoxy_chain_length = 0

        # For each fragment, count the number of carbon atoms
        for frag in fragments:
            carbon_atoms = [atom for atom in frag.GetAtoms() if atom.GetAtomicNum() == 6]
            c_count = len(carbon_atoms)
            # Determine if fragment is acyl or alkoxy based on attachment to dummy atom
            # Dummy atom (atomic number 0) replaces the bond we broke
            dummy_neighbors = [atom for atom in frag.GetAtoms() if atom.GetAtomicNum() == 0]
            if len(dummy_neighbors) != 1:
                continue  # Fragment should have one dummy atom
            dummy_atom = dummy_neighbors[0]
            connected_atoms = dummy_atom.GetNeighbors()
            if connected_atoms:
                neighbor = connected_atoms[0]
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() in [c_atom_idx, o_atom_idx]:
                    # This fragment contains the original carbonyl carbon
                    acyl_chain_length = c_count
                else:
                    alkoxy_chain_length = c_count

        # Check if both chains are long enough (14 carbons or more)
        if acyl_chain_length >= 14 and alkoxy_chain_length >= 14:
            return True, f"Ester with acyl chain length {acyl_chain_length} and alkoxy chain length {alkoxy_chain_length}"

    return False, "Ester present, but chains are not long enough for wax"

__metadata__ = {
    'chemical_class': {
        'name': 'wax',
        'definition': 'A chemical substance that is an organic compound or mixture of compounds that is composed of long-chain molecules and is malleable at ambient temperatures.'
    }
}
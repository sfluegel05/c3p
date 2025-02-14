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

    # Check for alcohol group (-OH)
    OH_pattern = Chem.MolFromSmarts('[OX2H]')
    OH_matches = mol.GetSubstructMatches(OH_pattern)
    if not OH_matches:
        return False, "No alcohol group (-OH) found"
    
    # Assume primary alcohol; get carbon attached to OH
    OH_atom_idx = OH_matches[0][0]
    OH_atom = mol.GetAtomWithIdx(OH_atom_idx)
    neighbors = OH_atom.GetNeighbors()
    if not neighbors:
        return False, "OH group is not connected to any carbon"
    # Carbon attached to OH group
    start_atom = neighbors[0]
    
    visited = set()
    chain_atoms = []

    # Function to recursively traverse the chain
    def traverse_chain(atom, prev_atom_idx):
        visited.add(atom.GetIdx())
        chain_atoms.append(atom)
        neighbors = [a for a in atom.GetNeighbors() if a.GetIdx() != prev_atom_idx]
        if len(neighbors) > 1:
            # Branch detected
            return False
        elif len(neighbors) == 0:
            # End of chain
            return True
        else:
            next_atom = neighbors[0]
            if next_atom.GetAtomicNum() != 6:
                # Non-carbon atom in chain
                return False
            return traverse_chain(next_atom, atom.GetIdx())
    
    # Start traversing from the carbon attached to OH
    is_linear = traverse_chain(start_atom, OH_atom_idx)
    if not is_linear:
        return False, "Chain is branched or contains non-carbon atoms"
    
    # Exclude the terminal OH carbon from chain length
    chain_length = len(chain_atoms)
    
    if chain_length < 1:
        return False, "Chain is too short"
    
    # Exclude the OH-attached carbon from the chain for counting isoprene units
    isoprene_chain = chain_atoms[1:]
    n_carbons = len(isoprene_chain)

    # Number of carbons should be a multiple of 5
    if n_carbons % 5 != 0:
        return False, f"Chain length ({n_carbons}) is not a multiple of 5"

    n_units = n_carbons // 5
    # Check each 5-carbon segment for isoprene pattern
    for i in range(n_units):
        segment_atoms = isoprene_chain[i*5:(i+1)*5]
        if len(segment_atoms) != 5:
            return False, "Incomplete isoprene unit detected"

        atom1 = segment_atoms[0]
        atom2 = segment_atoms[1]
        atom3 = segment_atoms[2]
        atom4 = segment_atoms[3]
        atom5 = segment_atoms[4]

        # Check that atom2 has a methyl substituent
        methyl_found = False
        for neighbor in atom2.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in [atom1.GetIdx(), atom3.GetIdx()]:
                # Check if it's a methyl group
                if len(neighbor.GetNeighbors()) == 1:
                    methyl_found = True
                    break
        if not methyl_found:
            return False, "Isoprene unit missing methyl group at correct position"
        
        # Check for double bonds between atom1-atom2 and/or atom3-atom4
        bond12 = mol.GetBondBetweenAtoms(atom1.GetIdx(), atom2.GetIdx())
        bond34 = mol.GetBondBetweenAtoms(atom3.GetIdx(), atom4.GetIdx())
        if bond12 is None or bond34 is None:
            return False, "Isoprene unit missing bonds"
        bond12_order = bond12.GetBondType()
        bond34_order = bond34.GetBondType()
        if bond12_order != Chem.rdchem.BondType.DOUBLE and bond34_order != Chem.rdchem.BondType.DOUBLE:
            return False, "Isoprene unit missing double bonds at correct positions"
        
    return True, "Molecule is a prenol with linear chain of isoprene units ending with -OH group"
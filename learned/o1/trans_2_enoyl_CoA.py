"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
"""
Classifies: trans-2-enoyl-CoA
"""
from rdkit import Chem

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a trans-2-enoyl-CoA based on its SMILES string.
    A trans-2-enoyl-CoA is an unsaturated fatty acyl-CoA that results from the formal condensation
    of the thiol group of coenzyme A with the carboxy group of any 2,3-trans-enoic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trans-2-enoyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define simplified CoA pattern focusing on key features
    coa_smarts = """
    [#7]-[#6]-[#6](=O)-[#6]-[#7]-[#6](=O)-[#6@H](-O)-[#6](-[#6])-[#6]-O-P(=O)(O)-O-P(=O)(O)-O
    """
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"

    # Define thioester linkage pattern
    thioester_smarts = "C(=O)S"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    # Get indices of CoA atoms
    coa_matches = mol.GetSubstructMatches(coa_pattern)
    coa_atom_indices = set(idx for match in coa_matches for idx in match)

    # For each thioester linkage, check if sulfur is connected to CoA
    for match in thioester_matches:
        carbonyl_c_idx, sulfur_idx = match[0], match[1]
        # Check if sulfur atom is part of CoA
        if sulfur_idx not in coa_atom_indices:
            continue  # Sulfur not connected to CoA

        # Extract acyl chain attached to carbonyl carbon
        acyl_chain_atoms = set()
        visited = set()
        stack = [carbonyl_c_idx]
        while stack:
            atom_idx = stack.pop()
            if atom_idx in visited or atom_idx in coa_atom_indices:
                continue
            visited.add(atom_idx)
            acyl_chain_atoms.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx in visited or neighbor_idx in coa_atom_indices:
                    continue
                stack.append(neighbor_idx)

        # Create sub-molecule of acyl chain
        acyl_chain = Chem.PathToSubmol(mol, acyl_chain_atoms)

        # Check that acyl chain has at least 3 carbons
        acyl_carbons = [atom for atom in acyl_chain.GetAtoms() if atom.GetAtomicNum() == 6]
        if len(acyl_carbons) < 3:
            continue  # Acyl chain too short

        # Assign stereochemistry to acyl chain
        Chem.AssignStereochemistry(acyl_chain, cleanIt=True, force=True)

        # Identify C1 (carbonyl carbon), C2, and C3 in acyl chain
        # Map original atom indices to acyl_chain atom indices
        atom_mapping = {}
        for atom in acyl_chain.GetAtoms():
            orig_idx = atom.GetProp('_ori_index')
            atom_mapping[int(orig_idx)] = atom.GetIdx()

        carbonyl_c_ac_idx = atom_mapping.get(carbonyl_c_idx, None)
        if carbonyl_c_ac_idx is None:
            continue  # Cannot find carbonyl carbon in acyl chain

        c1_atom = acyl_chain.GetAtomWithIdx(carbonyl_c_ac_idx)

        # Get C2
        neighbors_c1 = [nbr for nbr in c1_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if not neighbors_c1:
            continue  # No C2 found
        c2_atom = neighbors_c1[0]

        # Get C3
        neighbors_c2 = [nbr for nbr in c2_atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != c1_atom.GetIdx()]
        if not neighbors_c2:
            continue  # No C3 found
        c3_atom = neighbors_c2[0]

        # Check for double bond between C2 and C3
        bond_c2_c3 = acyl_chain.GetBondBetweenAtoms(c2_atom.GetIdx(), c3_atom.GetIdx())
        if bond_c2_c3 is None or bond_c2_c3.GetBondType() != Chem.rdchem.BondType.DOUBLE:
            continue  # No double bond between C2 and C3

        # Check if the double bond is trans (E configuration)
        stereo = bond_c2_c3.GetStereo()
        if stereo != Chem.rdchem.BondStereo.STEREOE:
            return False, "Double bond between C2 and C3 is not trans"

        # All conditions met
        return True, "Contains CoA moiety with acyl chain having trans double bond between C2 and C3"

    # If no suitable acyl chain found
    return False, "Acyl chain with required features not found"

__metadata__ = {
    'chemical_class': {
        'id': '',
        'name': 'trans-2-enoyl-CoA',
        'definition': 'An unsaturated fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any 2,3-trans-enoic acid.',
        'parents': []
    },
    'message': None,
    'success': True,
    'error': '',
}
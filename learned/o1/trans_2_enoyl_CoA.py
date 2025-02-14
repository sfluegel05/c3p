"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
"""
Classifies: trans-2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

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

    # Generate 3D coordinates to perceive stereochemistry
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

    # Define key substructures of CoA
    # Adenosine monophosphate part
    adenine_smarts = "n1cnc2c(ncnc12)"
    ribose_smarts = "OC[C@H]1O[C@H](n2cnc3c(ncnc32))[C@H](O)[C@@H]1O"
    phosphate_smarts = "OP(=O)(O)O"

    # Pantetheine arm with thiol group (we'll look for the thiol converted to thioester)
    pantetheine_smarts = "NCC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)O"

    # Combine to define simplified CoA pattern
    coa_pattern_smarts = f"{pantetheine_smarts}.{adenine_smarts}"

    coa_pattern = Chem.MolFromSmarts(coa_pattern_smarts)
    if coa_pattern is None:
        return False, "Failed to generate CoA substructure"

    # Check for CoA substructure
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"

    # Define thioester linkage pattern: carbonyl carbon connected to sulfur
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    # Find the acyl chain attached to the thioester
    for match in thioester_matches:
        carbonyl_c_idx = match[0]
        sulfur_idx = match[1]

        # Check that sulfur is part of CoA by ensuring it's in the CoA match
        coa_matches = mol.GetSubstructMatches(coa_pattern)
        if not coa_matches:
            continue  # Should not happen, CoA already checked
        coa_atom_indices = set(idx for match in coa_matches for idx in match)
        if sulfur_idx not in coa_atom_indices:
            continue  # Sulfur not part of CoA

        # Traverse the acyl chain starting from the carbonyl carbon
        acyl_chain_atoms = set()
        visited = set()
        stack = [carbonyl_c_idx]
        while stack:
            atom_idx = stack.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            acyl_chain_atoms.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx == sulfur_idx:
                    continue  # Skip sulfur atom (part of CoA)
                if neighbor_idx in coa_atom_indices:
                    continue  # Skip atoms in CoA
                if neighbor_idx not in visited:
                    stack.append(neighbor_idx)

        # Sort acyl chain atoms by connectivity starting from carbonyl carbon
        acyl_chain = Chem.PathToSubmol(mol, acyl_chain_atoms)
        acyl_chain = Chem.AddHs(acyl_chain)
        Chem.AssignStereochemistry(acyl_chain, cleanIt=True, force=True)

        # Check if acyl chain has at least 3 carbons
        acyl_carbons = [atom for atom in acyl_chain.GetAtoms() if atom.GetAtomicNum() == 6]
        if len(acyl_carbons) < 3:
            continue  # Not a fatty acyl chain

        # Get atoms for C1, C2, and C3
        carbonyl_c_ac = None
        c2_ac = None
        c3_ac = None

        # Find carbonyl carbon (C1)
        for atom in acyl_chain.GetAtoms():
            if atom.GetAtomicNum() == 6:
                for bond in atom.GetBonds():
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and bond.GetOtherAtom(atom).GetAtomicNum() == 8:
                        carbonyl_c_ac = atom
                        break
            if carbonyl_c_ac:
                break

        if carbonyl_c_ac is None:
            continue  # Cannot find carbonyl carbon

        # Traverse to find C2 and C3
        neighbors_c1 = [nbr for nbr in carbonyl_c_ac.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if not neighbors_c1:
            continue
        c2_ac = neighbors_c1[0]

        neighbors_c2 = [nbr for nbr in c2_ac.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != carbonyl_c_ac.GetIdx()]
        if not neighbors_c2:
            continue
        c3_ac = neighbors_c2[0]

        # Check for double bond between C2 and C3
        bond_c2_c3 = acyl_chain.GetBondBetweenAtoms(c2_ac.GetIdx(), c3_ac.GetIdx())
        if bond_c2_c3 is None or bond_c2_c3.GetBondType() != Chem.rdchem.BondType.DOUBLE:
            continue  # No double bond between C2 and C3

        # Check if the double bond between C2 and C3 is trans
        stereo = bond_c2_c3.GetStereo()
        if stereo != Chem.rdchem.BondStereo.STEREOE:
            return False, "Double bond between C2 and C3 is not trans"

        return True, "Contains CoA moiety with acyl chain having trans double bond between C2 and C3"

    return False, "Does not contain acyl chain with trans double bond between C2 and C3 attached to CoA"

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
"""
Classifies: CHEBI:87659 dodecanoate ester
"""
"""
Classifies: CHEBI:36000 dodecanoate ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester based on its SMILES string.
    A dodecanoate ester is an ester where the acyl chain (acid component) is derived from lauric acid 
    (dodecanoic acid, a saturated 12-carbon fatty acid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dodecanoate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify ester groups: pattern [C;D2](=O)[O;D1]
    ester_pattern = Chem.MolFromSmarts('[C;D2](=O)[O;D1][#6]')
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"

    # Iterate over each ester group
    for match in ester_matches:
        carbonyl_c_idx = match[0]  # Carbonyl carbon
        ester_o_idx = match[2]     # Ester oxygen

        # Break the bond between carbonyl carbon and ester oxygen
        bond_to_break = mol.GetBondBetweenAtoms(carbonyl_c_idx, ester_o_idx).GetIdx()
        fragmented_mol = Chem.FragmentOnBonds(mol, [bond_to_break], addDummies=False)

        # Get the fragments
        frags = Chem.GetMolFrags(fragmented_mol, asMols=True, sanitizeFrags=True)
        if not frags:
            continue  # No fragments obtained, skip to next match

        # The acyl fragment should contain the carbonyl carbon
        acyl_frag = None
        for frag in frags:
            atom_indices = [atom.GetIdx() for atom in frag.GetAtoms()]
            if any(atom.GetAtomicNum() == 6 and atom.GetIdx() == carbonyl_c_idx for atom in frag.GetAtoms()):
                acyl_frag = frag
                break
        if acyl_frag is None:
            continue  # Could not find acyl fragment, continue to next ester group

        # Analyze the acyl fragment
        num_carbons = 0
        num_double_bonds = 0
        has_rings = False
        is_linear = True

        for atom in acyl_frag.GetAtoms():
            if atom.GetAtomicNum() == 6:
                num_carbons += 1
                if atom.GetDegree() > 2 and not atom.IsInRing():
                    # Detect branching
                    is_linear = False
            for bond in atom.GetBonds():
                if bond.IsInRing():
                    has_rings = True
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and bond.GetIdx() != carbonyl_c_idx:
                    # Exclude the carbonyl double bond
                    num_double_bonds += 1

        if num_carbons != 12:
            continue  # Not a 12-carbon chain
        if has_rings:
            continue  # Contains rings
        if not is_linear:
            continue  # Contains branches
        if num_double_bonds != 1:
            continue  # Should have only one double bond (the carbonyl)

        # Passed all checks
        return True, "Contains dodecanoate ester group with lauric acid acyl chain"

    return False, "No dodecanoate ester groups with lauric acid found"
"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
"""
Classifies: trans-2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a trans-2-enoyl-CoA based on its SMILES string.
    A trans-2-enoyl-CoA is an unsaturated fatty acyl-CoA resulting from the formal condensation
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

    # Define Coenzyme A (CoA) moiety pattern
    coa_smarts = "[#7]-[#6]-[#6](=O)-[#6]-[#7]-[#6](=O)-[#6@H](-O)-[#6](-[#6])-[#6]-O-P(=O)(O)-O-P(=O)(O)-O"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if coa_pattern is None:
        return False, "Invalid CoA SMARTS pattern"

    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"

    # Define thioester linkage pattern (CoA-S-C(=O)-)
    thioester_smarts = "S-C(=O)-C"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    if thioester_pattern is None:
        return False, "Invalid thioester SMARTS pattern"

    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    # Check for trans double bond between C2 and C3 of acyl chain
    double_bond_found = False
    for thioester_match in thioester_matches:
        sulfur_idx, carbonyl_c_idx, c1_idx = thioester_match  # S, C(=O), C
        c1_atom = mol.GetAtomWithIdx(c1_idx)

        # Find adjacent carbon (C2)
        neighbors_c1 = [nbr for nbr in c1_atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != carbonyl_c_idx]
        if not neighbors_c1:
            continue  # No C2 found
        c2_atom = neighbors_c1[0]

        # Check for double bond between C2 and C3
        bonds = c2_atom.GetBonds()
        double_bond = None
        for bond in bonds:
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and bond.GetIdx() != c1_atom.GetIdx():
                double_bond = bond
                break
        if not double_bond:
            continue  # No double bond found

        # Check if the double bond is between C2 and C3
        c3_atom = double_bond.GetOtherAtom(c2_atom)

        # Check if double bond is trans (E)
        bond_stereo = double_bond.GetStereo()
        if bond_stereo != Chem.rdchem.BondStereo.STEREOE:
            continue  # Not a trans double bond

        double_bond_found = True
        break  # Suitable acyl chain found

    if not double_bond_found:
        return False, "Acyl chain does not have trans double bond between C2 and C3"

    return True, "Molecule is a trans-2-enoyl-CoA"

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
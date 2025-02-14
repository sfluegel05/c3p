"""
Classifies: CHEBI:17334 penicillin
"""
"""
Classifies: penicillin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_penicillin(smiles: str):
    """
    Determines if a molecule is a penicillin based on its SMILES string.
    A penicillin is defined as any substituted penam containing two methyl groups at position 2,
    a carboxylate group at position 3, and a carboxamido group at position 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a penicillin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define penam core with atom mapping
    penam_smarts = """
    [#7:1]-1-[#6:2]=[#7:3]-[#6:4](-[#6:5])(-[#6:6])-[#16:7]-1-[#6:8]-[#6:9]=[#8:10]
    """

    penam_core = Chem.MolFromSmarts(penam_smarts)
    if penam_core is None:
        return False, "Error in penam core SMARTS"

    match = mol.GetSubstructMatch(penam_core)
    if not match:
        return False, "Penam core not found"

    # Map atom indices from the match
    atom_indices = dict(zip([1,2,3,4,5,6,7,8,9,10], match))

    # Positions according to atom mapping
    pos2 = atom_indices[5]  # Carbon at position 2 (alpha to sulfur)
    pos3 = atom_indices[4]  # Carbon at position 3 (beta to sulfur)
    pos6 = atom_indices[2]  # Carbon at position 6 (attached to nitrogen)

    # Check for two methyl groups at position 2
    pos2_atom = mol.GetAtomWithIdx(pos2)
    methyl_count = 0
    for neighbor in pos2_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6 and neighbor.GetDegree() == 1:
            methyl_count += 1
    if methyl_count != 2:
        return False, f"Position 2 does not have two methyl groups, found {methyl_count}"

    # Check for carboxylate group at position 3
    pos3_atom = mol.GetAtomWithIdx(pos3)
    carboxylate_found = False
    for bond in pos3_atom.GetBonds():
        neighbor = bond.GetOtherAtom(pos3_atom)
        if neighbor.GetAtomicNum() == 6 and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            # Look for C(=O)[O-] or C(=O)O
            neighbor_idx = neighbor.GetIdx()
            substructure = Chem.MolFragmentToSmiles(mol, atomsToUse=[pos3, neighbor_idx], canonical=True)
            carboxylate_smarts = Chem.MolFromSmarts('C(=O)[O-,O]')
            if mol.HasSubstructMatch(carboxylate_smarts, useChirality=False):
                carboxylate_found = True
                break
    if not carboxylate_found:
        return False, "Carboxylate group at position 3 not found"

    # Check for carboxamido group at position 6
    pos6_atom = mol.GetAtomWithIdx(pos6)
    carboxamido_found = False
    for bond in pos6_atom.GetBonds():
        neighbor = bond.GetOtherAtom(pos6_atom)
        if neighbor.GetAtomicNum() == 7:
            # Look for NC(=O)
            amide_carbon = None
            for nb_bond in neighbor.GetBonds():
                nb_neighbor = nb_bond.GetOtherAtom(neighbor)
                if nb_neighbor.GetAtomicNum() == 6 and nb_bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    # Check if carbon has a =O
                    carbon = nb_neighbor
                    double_bonded_oxygen = False
                    for c_bond in carbon.GetBonds():
                        c_neighbor = c_bond.GetOtherAtom(carbon)
                        if c_neighbor.GetAtomicNum() == 8 and c_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            double_bonded_oxygen = True
                            break
                    if double_bonded_oxygen:
                        carboxamido_found = True
                        break
            if carboxamido_found:
                break
    if not carboxamido_found:
        return False, "Carboxamido group at position 6 not found"

    return True, "Molecule matches penicillin structure with required substituents"

__metadata__ = {
    'chemical_class': {
        'name': 'penicillin',
        'definition': 'Any member of the group of substituted penams containing two methyl substituents at position 2, a carboxylate substituent at position 3 and a carboxamido group at position 6.'
    }
}
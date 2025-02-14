"""
Classifies: CHEBI:26333 prostaglandin
"""
"""
Classifies: prostaglandin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    Prostaglandins are compounds derived from prostanoic acid (C20),
    containing specific functional groups and structural features.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check total number of carbon atoms (should be close to 20)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18 or c_count > 22:
        return False, f"Carbon count {c_count} not within prostaglandin range (18-22)"

    # Check for carboxylic acid group or its derivatives at one end
    carboxylic_group = Chem.MolFromSmarts('C(=O)[O;H1,-]')
    has_carboxylic = mol.HasSubstructMatch(carboxylic_group)
    ester_group = Chem.MolFromSmarts('C(=O)O[C]')
    has_ester = mol.HasSubstructMatch(ester_group)
    amide_group = Chem.MolFromSmarts('C(=O)N')
    has_amide = mol.HasSubstructMatch(amide_group)

    if not (has_carboxylic or has_ester or has_amide):
        return False, "No carboxylic acid group or its derivative found"

    # Check for hydroxyl groups
    hydroxyl_group = Chem.MolFromSmarts('[OX2H]')
    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_group)
    if not has_hydroxyl:
        return False, "No hydroxyl group found"

    # Ketone group is optional
    ketone_group = Chem.MolFromSmarts('C(=O)[#6]')
    has_ketone = mol.HasSubstructMatch(ketone_group)
    # Note: We won't enforce ketone presence

    # Check for cyclopentane ring (may or may not be present)
    ring_info = mol.GetRingInfo()
    has_cyclopentane = False
    for ring in ring_info.AtomRings():
        if len(ring) == 5:
            # Check if all atoms in the ring are carbons
            if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                has_cyclopentane = True
                break  # At least one cyclopentane ring found

    # Prostaglandins may have conjugated double bonds in the chain
    conjugated_diene = Chem.MolFromSmarts('C=C-C=C')
    has_conjugated_diene = mol.HasSubstructMatch(conjugated_diene)

    # Check for long aliphatic chain (at least 7 carbons in a row)
    aliphatic_chain = Chem.MolFromSmarts('[C;!R](-[C;!R]){6,}')
    has_long_chain = mol.HasSubstructMatch(aliphatic_chain)
    if not has_long_chain:
        return False, "No long aliphatic chain found"

    # Overall check: molecule should have either cyclopentane ring or conjugated diene
    if not (has_cyclopentane or has_conjugated_diene):
        return False, "Neither cyclopentane ring nor conjugated diene found"

    # If all checks pass
    return True, "Matches prostaglandin structural features"

# Examples of usage
if __name__ == "__main__":
    smiles_list = [
        # Examples of prostaglandins
        "CCCCCC[C@H](O)\\C=C\\C1=C(C\\C=C/CCCC(O)=O)C(=O)CC1",  # Prostaglandin B2
        "CCCCCCCCC(=O)CC[C@H]1[C@H](O)C[C@H](O)[C@@H]1C\\C=C/CCCC(=O)OC(C)C",  # Isopropyl unoprostone
        "C(CCC(O)=O)/C=C\\C[C@@H]1\\C(=C\\C=C/C=C\\CC)[C@@H](O)CC1=O",  # 15-deoxy-Δ12,14-prostaglandin J2
        # Non-prostaglandin for testing
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O"  # Nonacosanoic acid
    ]
    for smiles in smiles_list:
        result, reason = is_prostaglandin(smiles)
        print(f"SMILES: {smiles}\nIs prostaglandin: {result}\nReason: {reason}\n")
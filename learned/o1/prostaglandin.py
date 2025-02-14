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
    containing a cyclopentane ring with specific functional groups.

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
    carboxylic_group = Chem.MolFromSmarts('C(=O)[O-]')
    if carboxylic_group is None:
        return False, "Invalid carboxylic acid SMARTS pattern"
    has_carboxylic = mol.HasSubstructMatch(carboxylic_group)

    # Also consider carboxylic acid in protonated form
    carboxylic_group_protonated = Chem.MolFromSmarts('C(=O)O[H]')
    if carboxylic_group_protonated is None:
        return False, "Invalid protonated carboxylic acid SMARTS pattern"
    has_carboxylic_protonated = mol.HasSubstructMatch(carboxylic_group_protonated)

    # Check for ester group (e.g., in ester derivatives)
    ester_group = Chem.MolFromSmarts('C(=O)OC')
    if ester_group is None:
        return False, "Invalid ester SMARTS pattern"
    has_ester = mol.HasSubstructMatch(ester_group)

    if not (has_carboxylic or has_carboxylic_protonated or has_ester):
        return False, "No carboxylic acid group or ester derivative found"

    # Check for hydroxyl groups attached to the cyclopentane ring
    hydroxyl_group = Chem.MolFromSmarts('[C;R][C;R][C;R][C;R][C;R]')  # Cyclopentane ring
    if hydroxyl_group is None:
        return False, "Invalid cyclopentane ring SMARTS pattern"
    matches = mol.GetSubstructMatches(hydroxyl_group)
    has_cyclopentane = len(matches) > 0

    if not has_cyclopentane:
        return False, "Cyclopentane ring not found"

    # Check for hydroxyl groups on the ring
    ring_atoms = set()
    for match in matches:
        ring_atoms.update(match)
    has_hydroxyl_on_ring = False
    for atom_idx in ring_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetSymbol() == 'O' and nbr.GetDegree() == 1:
                has_hydroxyl_on_ring = True
                break
        if has_hydroxyl_on_ring:
            break

    if not has_hydroxyl_on_ring:
        return False, "No hydroxyl group attached to cyclopentane ring"

    # Check for double bonds in side chains
    double_bond = Chem.MolFromSmarts('C=C')
    if double_bond is None:
        return False, "Invalid double bond SMARTS pattern"
    has_double_bond = mol.HasSubstructMatch(double_bond)
    if not has_double_bond:
        return False, "No double bonds found in side chains"

    # If all checks pass
    return True, "Matches prostaglandin structural features"

# Examples of usage
if __name__ == "__main__":
    smiles_list = [
        # Examples of prostaglandins
        "CCCCC[C@H](O)\\C=C\\[C@H]1[C@H](O)C[C@H](O)[C@H]1C\\C=C/CC(O)=O",  # Prostaglandin F2alpha
        "CCCCCCCCC(=O)CC[C@H]1[C@H](O)C[C@H](O)[C@@H]1C\\C=C/CCCC(=O)OC(C)C",  # Isopropyl unoprostone
        "C(CCC(O)=O)/C=C\\C[C@@H]1\\C(=C\\C=C/C=C\\CC)[C@@H](O)CC1=O",  # 15-deoxy-Δ12,14-prostaglandin J2
        # Non-prostaglandin for testing
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O"  # Nonacosanoic acid
    ]
    for smiles in smiles_list:
        result, reason = is_prostaglandin(smiles)
        print(f"SMILES: {smiles}\nIs prostaglandin: {result}\nReason: {reason}\n")
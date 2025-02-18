"""
Classifies: CHEBI:26333 prostaglandin
"""
"""
Classifies: CHEBI:15846 prostaglandin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    Prostaglandins are derived from prostanoic acid, featuring a cyclopentane ring with two substituent chains,
    one ending in a carboxylic acid derivative, and containing oxygen functional groups and double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for cyclopentane ring (all carbons)
    cyclopentane = Chem.MolFromSmarts("C1CCCC1")
    if not mol.HasSubstructMatch(cyclopentane):
        return False, "No cyclopentane ring"

    # Check for at least two substituent chains on the cyclopentane
    matches = mol.GetSubstructMatches(cyclopentane)
    has_two_substituents = False
    for ring in matches:
        substituent_count = 0
        ring_atoms = set(ring)
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in ring_atoms:
                    substituent_count += 1
                    break  # Count each ring atom with at least one substituent
        if substituent_count >= 2:
            has_two_substituents = True
            break
    if not has_two_substituents:
        return False, "Cyclopentane ring does not have at least two substituent chains"

    # Check for carboxylic acid derivative (acid, ester, amide)
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[O,N]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid derivative group (acid, ester, amide)"

    # Check oxygen count (prostaglandins have multiple O-containing groups)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, f"Insufficient oxygen atoms ({o_count})"

    # Check carbon count (prostanoic acid is C20, derivatives may vary)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (15 <= c_count <= 25):
        return False, f"Carbon count {c_count} outside typical range (15-25)"

    # Check for at least one hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl groups present"

    # Check for at least one double bond
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "No double bonds present"

    return True, "Contains cyclopentane ring with two substituent chains, carboxylic acid derivative, hydroxyl groups, and double bonds"
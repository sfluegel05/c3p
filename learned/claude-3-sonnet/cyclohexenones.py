"""
Classifies: CHEBI:48953 cyclohexenones
"""
"""
Classifies: CHEBI:33738 cyclohexenone
Any six-membered alicyclic ketone having one double bond in the ring.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_cyclohexenone(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclohexenone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all 6-membered alicyclic rings
    ring_info = mol.GetRingInfo()
    alicyclic_6_rings = [ring for ring in ring_info.AtomRings() if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring)]

    if not alicyclic_6_rings:
        return False, "No 6-membered alicyclic rings found"

    # Check if any of the 6-membered alicyclic rings has exactly one C=O and one C=C in the ring
    for ring in alicyclic_6_rings:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        carbonyl_pattern = Chem.MolFromSmarts("[C=O]")
        carbonyl_matches = [mol.GetAtomWithIdx(match) for match in mol.GetSubstructMatches(carbonyl_pattern) if mol.GetAtomWithIdx(match) in ring_atoms]
        double_bond_pattern = Chem.MolFromSmarts("[C=C]")
        double_bond_matches = [mol.GetBondBetweenAtoms(match[0], match[1]) for match in mol.GetSubstructMatches(double_bond_pattern) if mol.GetBondBetweenAtoms(match[0], match[1]).IsInRing() and mol.GetAtomWithIdx(match[0]) in ring_atoms and mol.GetAtomWithIdx(match[1]) in ring_atoms]

        if len(carbonyl_matches) == 1 and len(double_bond_matches) == 1:
            # Check if C=O and C=C are in conjugation (separated by 2 bonds)
            carbonyl_idx = carbonyl_matches[0].GetIdx()
            double_bond_idx = [bond.GetBeginAtomIdx() for bond in double_bond_matches][0]
            path = Chem.GetShortestPath(mol.GetTopology(), carbonyl_idx, double_bond_idx)
            if len(path) == 3:
                return True, "Molecule is a cyclohexenone"

    return False, "No valid cyclohexenone rings found"
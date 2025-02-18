"""
Classifies: CHEBI:28963 amino sugar
"""
"""
Classifies: CHEBI:45685 amino sugar
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmiles

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    An amino sugar is a sugar with one or more alcoholic hydroxy groups replaced by amino groups (substituted or unsubstituted).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino sugar, False otherwise
        str: Reason for classification
    """
    mol = MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for presence of at least one amino group (any nitrogen atom)
    has_amino = any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms())
    if not has_amino:
        return False, "No amino groups found"

    # Check for sugar-like ring structure (5/6-membered with oxygen and multiple hydroxyls)
    rings = mol.GetRingInfo().AtomRings()
    sugar_ring_found = False
    for ring in rings:
        if len(ring) not in {5, 6}:
            continue
        # Check for oxygen in the ring
        if not any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring):
            continue
        # Count hydroxyl groups on ring carbons
        hydroxyl_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:  # Skip the ring oxygen
                continue
            # Check for hydroxyl groups (-OH) attached to ring carbons
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >= 1 and neighbor.GetDegree() == 1:
                    hydroxyl_count += 1
        if hydroxyl_count >= 2:
            sugar_ring_found = True
            break

    if not sugar_ring_found:
        return False, "No sugar ring with multiple hydroxyl groups found"

    return True, "Contains amino group and sugar ring structure with multiple hydroxyls"
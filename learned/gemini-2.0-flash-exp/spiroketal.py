"""
Classifies: CHEBI:72600 spiroketal
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    A spiroketal is a cyclic ketal in which the ketal carbon is the only common atom of two rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a spiroketal, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Identify ketal carbons: carbon bonded to two oxygens (sp3 hybridized), zero hydrogens and two carbons
    ketal_pattern = Chem.MolFromSmarts("[CX4H0](-[OX2])(-[OX2])([CX4])([CX4])")
    ketal_matches = mol.GetSubstructMatches(ketal_pattern)
    if not ketal_matches:
        return False, "No ketal carbon found."

    for match in ketal_matches:
        ketal_carbon_idx = match[0]
        ketal_carbon = mol.GetAtomWithIdx(ketal_carbon_idx)

         # Get oxygen atoms bonded to ketal carbon
        oxygen_indices = [neighbor.GetIdx() for neighbor in ketal_carbon.GetNeighbors() if neighbor.GetAtomicNum() == 8]
        if len(oxygen_indices) != 2:
            continue # Should not happen

        oxy1_idx = oxygen_indices[0]
        oxy2_idx = oxygen_indices[1]

        oxy1 = mol.GetAtomWithIdx(oxy1_idx)
        oxy2 = mol.GetAtomWithIdx(oxy2_idx)

        # 2. Verify that both oxygens are part of rings.
        if not oxy1.IsInRing() or not oxy2.IsInRing():
            continue

        # 3. Find all rings for each oxygen, excluding the ketal carbon.
        ring_info = mol.GetRingInfo()
        rings_oxy1 = [set(ring) for ring in ring_info.AtomRings() if oxy1_idx in ring]
        rings_oxy2 = [set(ring) for ring in ring_info.AtomRings() if oxy2_idx in ring]

        if not rings_oxy1 or not rings_oxy2:
             continue

        is_spiroketal_candidate = True
        #Check all ring combinations for intersection
        for ring1 in rings_oxy1:
            ring1.discard(ketal_carbon_idx)
            for ring2 in rings_oxy2:
                ring2.discard(ketal_carbon_idx)
                if ring1.intersection(ring2):
                    is_spiroketal_candidate = False
                    break
            if not is_spiroketal_candidate:
                break

        if is_spiroketal_candidate:
            return True, "Spiroketal structure identified."

    return False, "No spiroketal structure detected."
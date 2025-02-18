"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI:3-hydroxy fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_hydroxy_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA(4-) based on its SMILES string.
    The structure must contain a 3-hydroxy acyl thioester group, three phosphate groups with four deprotonated oxygens,
    and have a sufficiently long fatty acid chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches criteria, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for 3-hydroxy acyl thioester pattern (S-C(=O)-CH2-CH(OH)-*)
    thioester_pattern = Chem.MolFromSmarts("[S]-C(=O)-[CH2]-[CH](-O)")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Missing 3-hydroxy acyl thioester group"

    # Check acyl chain length (at least 8 carbons in total)
    # Starting from the carbonyl carbon (index 1 in the match), count the chain
    min_chain_length = 8
    for match in thioester_matches:
        s_idx = match[0]
        co_idx = match[1]
        c1_idx = match[2]
        c2_idx = match[3]  # CH(OH)

        # Traverse the chain from c2 onwards
        def count_carbons(atom_idx, visited):
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() != 6:
                return 0
            count = 1
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                    count += count_carbons(neighbor.GetIdx(), visited)
            return count

        total_carbons = count_carbons(c2_idx, set()) + 3  # include the first three carbons (C=O, CH2, CH(OH))
        if total_carbons >= min_chain_length:
            break
    else:
        return False, f"Acyl chain too short ({total_carbons} carbons), need at least {min_chain_length}"

    # Check for three phosphorus atoms (three phosphate groups)
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(p_atoms) != 3:
        return False, f"Found {len(p_atoms)} phosphorus atoms, need exactly 3"

    # Count deprotonated oxygens (formal charge -1) attached to phosphorus
    o_minus_count = 0
    for p in p_atoms:
        for neighbor in p.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetFormalCharge() == -1:
                o_minus_count += 1
    if o_minus_count != 4:
        return False, f"Found {o_minus_count} deprotonated phosphate oxygens, need exactly 4"

    # Verify overall charge state (may not be reliable if SMILES lacks explicit charges)
    total_charge = Chem.GetFormalCharge(mol)
    if total_charge != -4:
        return False, f"Formal charge is {total_charge}, should be -4"

    return True, "Contains 3-hydroxy acyl thioester, three phosphates with four deprotonated oxygens, and charge -4"
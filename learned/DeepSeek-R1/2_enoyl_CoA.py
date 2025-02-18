"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
"""
Classifies: 2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2-enoyl-CoA based on its SMILES string.
    A 2-enoyl-CoA has an unsaturated fatty acyl group attached to CoA via a thioester,
    with a double bond between positions 2 and 3 of the acyl chain (counting from the carbonyl).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-enoyl-CoA, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for CoA structure (pantetheine part)
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing CoA moiety"

    # Find thioester groups (C(=O)S-)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group"

    # Check acyl chain for double bond between positions 2 and 3 (from carbonyl)
    # SMARTS pattern: [CX3](=[OX1])[SX2][CX4H2][CX4H2][CX3]=[CX3]
    # This matches C(=O)S-C-C-C=C (double bond at position 3-4 from sulfur)
    acyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2][CX4H2][CX4H2][CX3]=[CX3]")
    if mol.HasSubstructMatch(acyl_pattern):
        return True, "Double bond between positions 2 and 3 of acyl chain"

    # Alternative check using atom traversal for conjugated/branched systems
    for (carbonyl_idx, sulfur_idx) in thioester_matches:
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Get the acyl chain atom bonded to carbonyl (excluding sulfur)
        acyl_start = None
        for neighbor in carbonyl_atom.GetNeighbors():
            if neighbor.GetIdx() != sulfur_idx and neighbor.GetSymbol() == 'C':
                acyl_start = neighbor
                break
        if not acyl_start:
            continue

        # Traverse two carbons along the acyl chain
        current = acyl_start
        position = 1
        next_carbon = None
        while position < 2:  # Move to position 2 (third carbon from carbonyl)
            next_carbons = [n for n in current.GetNeighbors() 
                            if n.GetSymbol() == 'C' and n.GetIdx() != carbonyl_idx]
            if not next_carbons:
                break
            next_carbon = next_carbons[0]
            current = next_carbon
            position += 1

        if position < 2:
            continue

        # Check for double bond between position 2 and 3
        for bond in current.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                other = bond.GetOtherAtom(current)
                if other.GetSymbol() == 'C' and other.GetIdx() != acyl_start.GetIdx():
                    return True, "Double bond between positions 2 and 3 of acyl chain"

    return False, "No double bond between positions 2 and 3 of acyl chain"
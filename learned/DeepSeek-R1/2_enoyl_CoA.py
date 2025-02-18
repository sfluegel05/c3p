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
    # SMARTS pattern: [CX3](=[OX1])[SX2][CX4][CX3]=[CX3]
    # This matches C(=O)S-C-C=C (double bond at positions 2-3 of acyl chain)
    acyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2][CX4][CX3]=[CX3]")
    if mol.HasSubstructMatch(acyl_pattern):
        return True, "Double bond between positions 2 and 3 of acyl chain"

    # Alternative check using atom traversal for conjugated/branched systems
    for match in thioester_matches:
        # match contains (carbonyl_c, O, S)
        sulfur_idx = match[2]
        sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)
        
        # Find the acyl chain carbon attached to sulfur (excluding the carbonyl)
        acyl_start = None
        for neighbor in sulfur_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != match[0]:  # not the carbonyl carbon
                acyl_start = neighbor
                break
        if not acyl_start:
            continue
        
        # Traverse two carbons along the acyl chain from sulfur
        current = acyl_start
        position = 1
        found_double = False
        
        # Check next carbon (position 2) for double bond to position 3
        next_c = None
        for neighbor in current.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != sulfur_idx:
                next_c = neighbor
                break
        if not next_c:
            continue
        
        # Check bonds from next_c for double bond
        for bond in next_c.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetOtherAtomIdx(next_c.GetIdx()) != current.GetIdx():
                # Found double bond between positions 2 and 3
                return True, "Double bond between positions 2 and 3 of acyl chain"

    return False, "No double bond between positions 2 and 3 of acyl chain"
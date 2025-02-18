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
    with a double bond between positions 2 and 3 of the acyl chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-enoyl-CoA, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for CoA structure (pantetheine part: SCCNC(=O)CCNC(=O)...)
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing CoA moiety"

    # Check for thioester group (C(=O)S)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group found"

    # Check acyl chain for double bond between positions 2 and 3
    # SMARTS: C(=O)S-C-C=*
    acyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2][CX4H2][CX3]=[CX3]")
    if mol.HasSubstructMatch(acyl_pattern):
        return True, "Acyl chain has double bond between positions 2 and 3"
    
    # Check if the double bond is in a different position but still between 2 and 3
    # (handling possible different configurations or branches)
    # Get the carbonyl carbon from thioester
    for match in thioester_matches:
        carbonyl_carbon = match[0]
        # Get the attached acyl chain carbon (next to carbonyl)
        neighbors = [n for n in mol.GetAtomWithIdx(carbonyl_carbon).GetNeighbors() 
                     if n.GetSymbol() != 'O' and n.GetSymbol() != 'S']
        if not neighbors:
            continue
        acyl_start = neighbors[0]
        # Traverse the acyl chain to find double bond between positions 2 and 3
        # Position 2 is the first carbon after acyl_start (position 1)
        # Position 3 is the next, so we need a double bond between them
        for bond in acyl_start.GetBonds():
            if bond.GetBeginAtomIdx() == acyl_start.GetIdx():
                next_atom = bond.GetEndAtom()
            else:
                next_atom = bond.GetBeginAtom()
            if bond.GetBondType() == Chem.BondType.DOUBLE and next_atom.GetSymbol() == 'C':
                return True, "Double bond found between positions 2 and 3 of acyl chain"
            # Check next bond in the chain
            for next_bond in next_atom.GetBonds():
                if next_bond.GetBondType() == Chem.BondType.DOUBLE and next_bond.GetOtherAtomIdx(next_atom.GetIdx()) != acyl_start.GetIdx():
                    return True, "Double bond found between positions 2 and 3 of acyl chain"

    return False, "No double bond between positions 2 and 3 of acyl chain"
"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
"""
Classifies: CHEBI:28645 unsaturated fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acyl-CoA based on its SMILES string.
    An unsaturated fatty acyl-CoA is characterized by a Coenzyme A moiety linked via a thioester
    bond to an unsaturated fatty acid (containing at least one double bond in the acyl chain).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for presence of thioester group (C(=O)S)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group detected"

    # Check for CoA structural motif (S-C-C-N-C=O pantetheine part)
    coa_pattern = Chem.MolFromSmarts("[SX2][CX4][CX4][NX3][CX3](=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"

    # Check each thioester group for unsaturated acyl chain
    for match in thioester_matches:
        carbonyl_idx = match[0]
        sulfur_idx = match[2]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)

        # Find acyl chain starting atom (non-O/S neighbor of carbonyl)
        acyl_start = None
        for neighbor in carbonyl_atom.GetNeighbors():
            if neighbor.GetAtomicNum() not in {8, 16}:  # Exclude O and S
                acyl_start = neighbor
                break
        if not acyl_start:
            continue  # Invalid thioester configuration

        # Traverse acyl chain to find any double bonds
        visited = {carbonyl_idx, sulfur_idx}
        stack = [acyl_start]
        acyl_atoms = set()
        while stack:
            atom = stack.pop()
            if atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())
            acyl_atoms.add(atom.GetIdx())
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in visited:
                    stack.append(neighbor)

        # Check for double bonds within the acyl chain
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                a1 = bond.GetBeginAtomIdx()
                a2 = bond.GetEndAtomIdx()
                if (a1 in acyl_atoms) and (a2 in acyl_atoms):
                    return True, "Contains unsaturated acyl chain attached to CoA via thioester"

    return False, "No unsaturated acyl chain found in thioester groups"
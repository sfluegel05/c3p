"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA based on its SMILES string.
    A monounsaturated fatty acyl-CoA has a CoA group attached via a thioester bond to a fatty acid
    with exactly one carbon-carbon double bond in the acyl chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for CoA thioester structure (C=O connected to S and CoA backbone)
    coa_thioester = Chem.MolFromSmarts('[CX3](=[OX1])SCCNC(=O)')
    matches = mol.GetSubstructMatches(coa_thioester)
    if not matches:
        return False, "Missing CoA thioester group"

    # Check each potential acyl chain (typically only one)
    for match in matches:
        carbonyl_c = mol.GetAtomWithIdx(match[0])
        sulfur = mol.GetAtomWithIdx(match[2])

        # Find acyl chain start (connected to carbonyl C, not O or S)
        acyl_start = None
        for neighbor in carbonyl_c.GetNeighbors():
            if neighbor.GetIdx() not in {match[1], sulfur.GetIdx()}:
                acyl_start = neighbor
                break
        if not acyl_start:
            continue

        # Traverse acyl chain to count C=C bonds
        visited = set([carbonyl_c.GetIdx(), sulfur.GetIdx()])
        stack = [acyl_start]
        visited.add(acyl_start.GetIdx())
        double_bonds = 0

        while stack:
            atom = stack.pop()
            for bond in atom.GetBonds():
                other = bond.GetOtherAtom(atom)
                if other.GetIdx() in visited:
                    continue
                # Count C=C bonds
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    if bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6:
                        double_bonds += 1
                visited.add(other.GetIdx())
                stack.append(other)

        if double_bonds == 1:
            return True, "Exactly one C=C bond in acyl chain"
        else:
            return False, f"Found {double_bonds} C=C bonds in acyl chain"

    return False, "No valid acyl chain with one double bond found"
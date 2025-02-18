"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA.
    Criteria: Contains CoA thioester, exactly one non-aromatic C=C bond in an acyclic aliphatic acyl chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if criteria met
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Match CoA thioester pattern (C=O connected to S and CoA backbone)
    thioester_pattern = Chem.MolFromSmarts('[CX3](=[OX1])SCCNC(=O)')
    matches = mol.GetSubstructMatches(thioester_pattern)
    if not matches:
        return False, "Missing CoA thioester"

    for match in matches:
        carbonyl_c = mol.GetAtomWithIdx(match[0])
        sulfur = mol.GetAtomWithIdx(match[2])

        # Find acyl chain starting carbon (connected to carbonyl, not O/S)
        acyl_start = None
        for neighbor in carbonyl_c.GetNeighbors():
            if neighbor.GetIdx() not in {match[1], sulfur.GetIdx()}:
                acyl_start = neighbor
                break
        if not acyl_start:
            continue

        # Traverse acyl chain to check structure
        visited = set([carbonyl_c.GetIdx(), sulfur.GetIdx()])
        stack = [acyl_start]
        visited.add(acyl_start.GetIdx())
        double_bonds = 0
        has_aromatic = False
        has_ring = False

        while stack and not (has_aromatic or has_ring):
            atom = stack.pop()
            
            # Check if atom is in a ring
            if atom.IsInRing():
                has_ring = True
                break
            
            for bond in atom.GetBonds():
                other = bond.GetOtherAtom(atom)
                if other.GetIdx() in visited:
                    continue
                
                # Check for aromatic bonds
                if bond.GetBondType() == Chem.BondType.AROMATIC:
                    has_aromatic = True
                    break
                
                # Count non-aromatic C=C bonds
                if (bond.GetBondType() == Chem.BondType.DOUBLE and
                    bond.GetBeginAtom().GetAtomicNum() == 6 and
                    bond.GetEndAtom().GetAtomicNum() == 6):
                    double_bonds += 1
                
                visited.add(other.GetIdx())
                stack.append(other)
        
        if has_aromatic or has_ring:
            continue  # Skip invalid acyl chains
        
        if double_bonds == 1:
            return True, "Exactly one non-aromatic C=C in acyclic aliphatic chain"
        else:
            return False, f"{double_bonds} C=C bonds in acyl chain"

    return False, "No valid acyl chain with one non-aromatic C=C"
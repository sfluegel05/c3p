"""
Classifies: CHEBI:64583 sphingomyelin
"""
"""
Classifies: CHEBI:17892 sphingomyelin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    Sphingomyelins have a sphingoid base with an amide-linked fatty acid and a phosphorylcholine group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingomyelin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphorylcholine group (O-P(=O)(OCC[N+](C)(C)C)-O- connected to carbon)
    pcho_pattern = Chem.MolFromSmarts("[O-]P(=O)(OCC[N+](C)(C)C)OC")
    pcho_matches = mol.GetSubstructMatches(pcho_pattern)
    if len(pcho_matches) != 1:
        return False, f"Found {len(pcho_matches)} phosphorylcholine groups, need exactly 1"

    # Check for amide-linked fatty acid (N-C(=O))
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) != 1:
        return False, f"Found {len(amide_matches)} amide groups, need exactly 1"

    # Verify amide is not part of the phosphorylcholine group
    amide_nitrogen_idx = amide_matches[0][0]
    pcho_atoms = set(pcho_matches[0])
    if amide_nitrogen_idx in pcho_atoms:
        return False, "Amide group is part of phosphorylcholine"

    # Check fatty acid chain length (minimum 12 carbons)
    carbonyl_carbon_idx = amide_matches[0][1]
    fatty_acid_carbons = 0
    # Traverse from carbonyl carbon, excluding the amide part
    visited = set()
    stack = [(carbonyl_carbon_idx, 0)]
    while stack:
        idx, depth = stack.pop()
        if idx in visited:
            continue
        visited.add(idx)
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6:
            continue
        fatty_acid_carbons = max(fatty_acid_carbons, depth + 1)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() == amide_nitrogen_idx:
                continue  # Skip back to amide nitrogen
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                stack.append((neighbor.GetIdx(), depth + 1))
    if fatty_acid_carbons < 12:
        return False, f"Fatty acid chain too short ({fatty_acid_carbons} carbons)"

    # Check sphingoid base structure (long chain with hydroxyl)
    sphingoid_carbon_idx = pcho_matches[0][-1]
    sphingoid_base_oxygen = pcho_matches[0][3]  # Oxygen connected to carbon
    # Check for adjacent hydroxyl group (part of sphingoid base)
    sphingoid_carbon = mol.GetAtomWithIdx(sphingoid_carbon_idx)
    for neighbor in sphingoid_carbon.GetNeighbors():
        if neighbor.GetAtomicNum() == 8 and neighbor.GetIdx() != sphingoid_base_oxygen:
            break
    else:
        # No additional hydroxyl found; check for double bond in sphingoid base
        # Look for double bonds in the sphingoid chain
        dbl_bonds = sum(1 for bond in sphingoid_carbon.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
        if dbl_bonds == 0:
            # Check if there's a double bond elsewhere in the sphingoid chain
            # Traverse the chain from sphingoid carbon
            visited = set()
            stack = [sphingoid_carbon_idx]
            has_double_bond = False
            while stack:
                idx = stack.pop()
                if idx in visited:
                    continue
                visited.add(idx)
                atom = mol.GetAtomWithIdx(idx)
                for bond in atom.GetBonds():
                    if bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetBeginAtomIdx() != sphingoid_base_oxygen and bond.GetEndAtomIdx() != sphingoid_base_oxygen:
                        has_double_bond = True
                        break
                    neighbor_idx = bond.GetOtherAtomIdx(idx)
                    if mol.GetAtomWithIdx(neighbor_idx).GetAtomicNum() == 6 and neighbor_idx not in visited:
                        stack.append(neighbor_idx)
                if has_double_bond:
                    break
            if not has_double_bond:
                return False, "Sphingoid base lacks hydroxyl or double bond"

    # Check molecular weight (typical sphingomyelins are >600 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da)"

    return True, "Contains phosphorylcholine, amide-linked fatty acid, and sphingoid base structure"
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

    # Improved phosphorylcholine pattern (more flexible atom order)
    pcho_pattern = Chem.MolFromSmarts("[O-]P(=O)([OX2])OCC[N+](C)(C)C")
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

    # Check fatty acid chain is hydrocarbon with minimum 12 carbons
    carbonyl_carbon = mol.GetAtomWithIdx(amide_matches[0][1])
    chain_carbons = 0
    visited = set()
    stack = [(carbonyl_carbon.GetIdx(), 0)]
    while stack:
        idx, depth = stack.pop()
        if idx in visited:
            continue
        visited.add(idx)
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6 or atom.GetFormalCharge() != 0:
            continue  # Must be neutral carbon
        chain_carbons = max(chain_carbons, depth + 1)
        for bond in atom.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetBeginAtomIdx() == idx and bond.GetEndAtom().GetAtomicNum() == 8:
                continue  # Skip ester bonds
            neighbor = bond.GetOtherAtomIdx(idx)
            if neighbor not in visited and mol.GetAtomWithIdx(neighbor).GetAtomicNum() == 6:
                stack.append((neighbor, depth + 1))
    if chain_carbons < 12:
        return False, f"Fatty acid chain too short ({chain_carbons} carbons)"

    # Check sphingoid base structure (long chain with at least one hydroxyl)
    # Find the sphingoid base connected to the phosphorylcholine oxygen
    pcho_o = [atom for atom in pcho_matches[0] if mol.GetAtomWithIdx(atom).GetAtomicNum() == 8 and mol.GetAtomWithIdx(atom).GetDegree() == 2][0]
    sphingoid_carbon = [n for n in mol.GetAtomWithIdx(pcho_o).GetNeighbors() if n.GetAtomicNum() == 6][0]
    
    # Traverse sphingoid base to find hydroxyl group and chain length
    sphingoid_oh = False
    visited = set()
    stack = [sphingoid_carbon.GetIdx()]
    while stack:
        idx = stack.pop()
        if idx in visited:
            continue
        visited.add(idx)
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1:
            sphingoid_oh = True
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                stack.append(neighbor.GetIdx())
    if not sphingoid_oh:
        return False, "Sphingoid base lacks hydroxyl group"

    return True, "Contains phosphorylcholine, amide-linked fatty acid, and sphingoid base structure"
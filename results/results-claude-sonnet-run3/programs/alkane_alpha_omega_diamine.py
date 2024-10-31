from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alkane_alpha_omega_diamine(smiles: str):
    """
    Determines if a molecule is an alkane-alpha,omega-diamine.
    These are primary diamines of the form H2NCH2(CH2)nCH2NH2 where n >= 0.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkane-alpha,omega-diamine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for exactly 2 nitrogen atoms
    nitrogens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if len(nitrogens) != 2:
        return False, "Must have exactly 2 nitrogen atoms"

    # Check that both nitrogens are primary amines (-NH2)
    for n in nitrogens:
        if n.GetTotalDegree() != 3:  # N should have 3 bonds total (2 H + 1 C)
            return False, "All nitrogen atoms must be primary amines (-NH2)"
        if len([x for x in n.GetNeighbors() if x.GetAtomicNum() == 6]) != 1:
            return False, "Each nitrogen must be bonded to exactly one carbon"

    # Get the shortest path between the two nitrogens
    path = Chem.GetShortestPath(mol, nitrogens[0].GetIdx(), nitrogens[1].GetIdx())
    if not path:
        return False, "No path found between nitrogen atoms"

    # Check that all atoms in the path (except N) are sp3 carbons
    for atom_idx in path[1:-1]:  # Exclude the terminal N atoms
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() != 6:
            return False, "Path between nitrogens must contain only carbon atoms"
        if atom.GetHybridization() != Chem.HybridizationType.SP3:
            return False, "All carbons must be sp3 hybridized"

    # Check that the carbons in the path only connect to H, OH or other carbons in path
    for atom_idx in path[1:-1]:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() not in [1, 6, 7, 8]:
                return False, "Carbon chain can only contain C, H, N, O atoms"
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in path:
                return False, "No branching allowed in carbon chain"
            if neighbor.GetAtomicNum() == 8:
                if not any(n.GetAtomicNum() == 1 for n in neighbor.GetNeighbors()):
                    return False, "Oxygen must be part of OH group"

    chain_length = len(path) - 2  # Subtract the two N atoms
    if chain_length < 1:
        return False, "Carbon chain must have at least one carbon"

    return True, f"Valid alkane-alpha,omega-diamine with {chain_length} carbons in chain"
# Pr=1.0
# Recall=0.8888888888888888
"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: Epoxy Fatty Acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    An epoxy fatty acid is a heterocyclic fatty acid containing an epoxide ring as part of its structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    ca_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not ca_matches:
        return False, "No carboxylic acid group found"
    # Get the index of the carbonyl carbon in the carboxylic acid group
    ca_c_index = ca_matches[0][0]

    # Check that the molecule has exactly one ring
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return False, "Molecule does not have exactly one ring (epoxide ring)"

    # Get ring atom indices and verify it's a three-membered ring
    ring_atoms = ring_info.AtomRings()[0]
    if len(ring_atoms) != 3:
        return False, "Ring is not a three-membered ring (epoxide ring)"

    # Check that the ring contains one oxygen and two carbons
    num_oxygen = 0
    num_carbon = 0
    for idx in ring_atoms:
        atom = mol.GetAtomWithIdx(idx)
        atomic_num = atom.GetAtomicNum()
        if atomic_num == 6:
            num_carbon += 1
        elif atomic_num == 8:
            num_oxygen += 1
        else:
            return False, "Ring contains atoms other than carbon and oxygen"
    if num_carbon != 2 or num_oxygen != 1:
        return False, "Ring is not an epoxide ring with two carbons and one oxygen"

    # Get indices of the epoxide carbons
    epoxide_c_indices = [idx for idx in ring_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]

    # Check for a path from the carboxylic acid carbon to the epoxide carbons
    # Ensure the path consists only of carbon atoms (typical of fatty acid chains)
    for c_idx in epoxide_c_indices:
        path = Chem.rdmolops.GetShortestPath(mol, ca_c_index, c_idx)
        if len(path) < 8:
            return False, "Chain between carboxylic acid and epoxide ring is too short"
        # Exclude the start and end atoms (already verified)
        path_atoms = [mol.GetAtomWithIdx(i) for i in path[1:-1]]
        for atom in path_atoms:
            if atom.GetAtomicNum() != 6:
                return False, "Path from carboxylic acid to epoxide contains heteroatoms"
            if atom.IsInRing():
                return False, "Path from carboxylic acid to epoxide contains rings"
            # Optionally, check for branching (atoms with more than 2 heavy atom neighbors)
            if len([nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]) > 2:
                return False, "Aliphatic chain is branched"
    # Ensure molecule contains only acceptable elements (C, H, O)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6, 8):
            return False, f"Unacceptable atom type found: atomic number {atom.GetAtomicNum()}"

    return True, "Molecule meets the criteria for an epoxy fatty acid"
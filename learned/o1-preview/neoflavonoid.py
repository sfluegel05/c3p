"""
Classifies: CHEBI:71971 neoflavonoid
"""
"""
Classifies: neoflavonoid
"""
from rdkit import Chem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.
    A neoflavonoid is any 1-benzopyran with an aryl substituent at position 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a neoflavonoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for 1-benzopyran core with atom mapping at position 4
    benzopyran_smarts = """
    [O]1[C;H2][C;H]([#6:4])[C;H]=[C;H][C;H]=[C;H]1
    """
    benzopyran_core = Chem.MolFromSmarts(benzopyran_smarts)
    if benzopyran_core is None:
        return False, "Invalid benzopyran SMARTS pattern"

    # Search for the benzopyran core in the molecule
    matches = mol.GetSubstructMatches(benzopyran_core)

    if not matches:
        return False, "1-benzopyran core not found"

    # For each match, check if there is an aryl substituent at position 4
    for match in matches:
        # Get the atom index for position 4 from the SMARTS mapping
        # In the SMARTS pattern, atom with map number 4 corresponds to position 4
        atom_idx_4 = match[benzopyran_core.GetSubstructMatch(Chem.MolFromSmarts("[#6:4]"))[0]]
        atom_4 = mol.GetAtomWithIdx(atom_idx_4)

        # Get neighbors of atom at position 4 that are not part of the benzopyran core
        benzopyran_atom_indices = set(match)
        neighbor_atoms = [a for a in atom_4.GetNeighbors() if a.GetIdx() not in benzopyran_atom_indices]

        # Check if any neighbor is part of an aromatic ring (aryl group)
        for neighbor in neighbor_atoms:
            if neighbor.IsAromatic():
                return True, "Molecule is a neoflavonoid with aryl substituent at position 4"

            # Check if the neighbor is connected to an aromatic system
            rings = neighbor.GetOwningMol().GetRingInfo().AtomRings()
            for ring in rings:
                if neighbor.GetIdx() in ring:
                    if all([mol.GetAtomWithIdx(i).IsAromatic() for i in ring]):
                        return True, "Molecule is a neoflavonoid with aryl substituent at position 4"

    return False, "No aryl substituent at position 4 of the 1-benzopyran core"
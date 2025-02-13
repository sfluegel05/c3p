"""
Classifies: CHEBI:16385 organic sulfide
"""
"""
Classifies: CHEBI:34842 organic sulfide
Organic sulfides are compounds with the structure RSR (R =/= H).
They were once called thioethers.
"""

from rdkit import Chem

def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organic sulfide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule contains any sulfur atoms
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]
    if not sulfur_atoms:
        return False, "No sulfur atoms found"

    # Check if at least one sulfur atom has exactly two non-hydrogen substituents
    has_organic_sulfide = False
    for s_atom in sulfur_atoms:
        non_h_neighbors = [nb for nb in s_atom.GetNeighbors() if nb.GetSymbol() != "H"]
        if len(non_h_neighbors) == 2:
            has_organic_sulfide = True
            break

    if has_organic_sulfide:
        return True, "Contains at least one sulfur atom with two non-hydrogen substituents (RSR)"
    else:
        return False, "No sulfur atoms with exactly two non-hydrogen substituents found"
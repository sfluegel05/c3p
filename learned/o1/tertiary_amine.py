"""
Classifies: CHEBI:32876 tertiary amine
"""
"""
Classifies: CHEBI:32876 tertiary amine
"""

from rdkit import Chem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.
    A tertiary amine has a nitrogen atom bonded to three carbon atoms and no hydrogens attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Flag to indicate if a tertiary amine nitrogen is found
    is_tertiary_amine = False

    # Iterate over all atoms in the molecule
    for atom in mol.GetAtoms():
        # Check if the atom is nitrogen
        if atom.GetAtomicNum() == 7:
            # Get the number of bonded atoms (degree)
            degree = atom.GetDegree()
            # Get the total number of hydrogens (implicit and explicit)
            total_h = atom.GetTotalNumHs()
            # Get formal charge
            formal_charge = atom.GetFormalCharge()

            # Check if nitrogen has degree 3, zero hydrogens, and no formal charge
            if degree == 3 and total_h == 0 and formal_charge == 0:
                # Check if all neighbors are carbon atoms
                neighbor_atoms = [nbr.GetAtomicNum() for nbr in atom.GetNeighbors()]
                if all(neighbor == 6 for neighbor in neighbor_atoms):
                    is_tertiary_amine = True
                    break  # No need to check further

    if is_tertiary_amine:
        return True, "Contains a nitrogen atom bonded to three carbon atoms with no hydrogens attached"
    else:
        return False, "No tertiary amine group found"
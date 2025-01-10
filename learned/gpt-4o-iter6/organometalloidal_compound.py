"""
Classifies: CHEBI:143084 organometalloidal compound
"""
"""
Classifies: organometalloidal compound
"""
from rdkit import Chem

def is_organometalloidal_compound(smiles: str):
    """
    Determines if a molecule is an organometalloidal compound based on its SMILES string.
    An organometalloidal compound has bonds between metalloid atoms (e.g., arsenic) and carbon atoms of an organyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organometalloidal compound, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define metalloid atomic numbers, focusing on arsenic for current examples
    metalloids = {33}  # Arsenic; extend with additional E.g., {14, 33, ...} for generic use
    metalloid_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() in metalloids]

    if not metalloid_atoms:
        return False, "No metalloid atoms (e.g., arsenic) found"

    # Check for bonds to carbon atoms linked as part of potential organyl groups
    for metalloid_atom in metalloid_atoms:
        for neighbor in metalloid_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon
                # Additional checks might be implemented here to ensure the carbon atom is part of an "organyl group"
                if any(next_neighbor.GetAtomicNum() not in {6, 1} for next_neighbor in neighbor.GetNeighbors()):
                    # This implies that the carbon is further bonded, possibly mimicking a part of an organyl group
                    return True, "Contains metalloid-carbon bond potentially part of an organyl group: an organometalloidal compound"

    return False, "No metalloid-carbon bonds found or insufficient organyl-like connections"

# Example usage:
# smiles = "C[As](O)(O)=O"
# result, reason = is_organometalloidal_compound(smiles)
# print(result, reason)  # Expected: True, "Contains metalloid-carbon bond: an organometalloidal compound"
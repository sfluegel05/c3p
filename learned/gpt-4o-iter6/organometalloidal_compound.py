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
    An organometalloidal compound has bonds between metalloid atoms (e.g., arsenic) and carbon atoms.

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

    # Define metalloid atomic numbers (arsenic, silicone, etc.)
    metalloids = {33}  # Add others if needed in future
    metalloid_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() in metalloids]

    if not metalloid_atoms:
        return False, "No metalloid atoms (e.g., arsenic) found"

    # Check for bonds to carbon atoms linked as part of organyl groups
    for metalloid_atom in metalloid_atoms:
        for neighbor in metalloid_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon
                # Further checks could be added here to validate the organyl nature
                # For now, we assume direct C bonding suffices as identified examples match
                return True, "Contains metalloid-carbon bond: an organometalloidal compound"

    return False, "No metalloid-carbon bonds found"

# Example usage:
# smiles = "C[As](O)(O)=O"
# result, reason = is_organometalloidal_compound(smiles)
# print(result, reason)  # Expected: True, "Contains metalloid-carbon bond: an organometalloidal compound"
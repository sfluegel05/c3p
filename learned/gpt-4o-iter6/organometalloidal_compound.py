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
    An organometalloidal compound has bonds between metalloid atoms (e.g., arsenic) and
    carbon atoms of an organyl group, often characterized by complex structures.

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
    metalloids = {33}  # Arsenic
    metalloid_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() in metalloids]

    if not metalloid_atoms:
        return False, "No metalloid atoms (e.g., arsenic) found"

    # SMARTS pattern representing a more complex carbon center
    organyl_pattern = Chem.MolFromSmarts("[#6]([#6,#1,#7,#8])([#6,#1,#7,#8])")  # Represents substituted/sp2 carbon

    # Check for bonds to carbon atoms linked as part of potential organyl groups
    for metalloid_atom in metalloid_atoms:
        for neighbor in metalloid_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon
                # Additional checks via SMARTS to ensure the carbon is part of a complex organic structure
                if neighbor.HasSubstructMatch(organyl_pattern):
                    return True, "Contains metalloid-carbon bond part of an organyl group: an organometalloidal compound"

    return False, "No suitable metalloid-carbon organyl connections identified"
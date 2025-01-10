"""
Classifies: CHEBI:25106 macrolide
"""
from rdkit import Chem

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.
    A macrolide is characterized by a macrocyclic lactone with a ring of twelve or more members derived from a polyketide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macrolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure the ring information is initialized
    ring_info = mol.GetRingInfo()
    if not ring_info or not ring_info.AtomRings():
        return False, "No rings identified in the molecule"

    # Define a SMARTS pattern that represents a basic form of a macrocyclic lactone
    lactone_pattern = Chem.MolFromSmarts('C1OC(=O)C1')

    # Iterate over all atom rings
    for ring_atoms in ring_info.AtomRings():
        # Check if the ring has 12 or more members
        if len(ring_atoms) >= 12:
            # Create a sub-molecule for the ring to examine for the lactone
            submol = Chem.PathToSubmol(mol, ring_atoms)
            # Check if this sub-molecule contains a lactone pattern
            if submol.HasSubstructMatch(lactone_pattern):
                return True, "Contains a macrocyclic lactone with 12 or more members"
    
    return False, "Does not contain a macrocyclic lactone with 12 or more members"
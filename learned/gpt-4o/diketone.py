"""
Classifies: CHEBI:46640 diketone
"""
from rdkit import Chem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone based on its SMILES string.
    A diketone contains exactly two ketone functionalities.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a diketone, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a refined SMARTS pattern for ketone functionality:
    # carbonyl group (C=O) bonded only to carbons, allowing no heteroatom
    ketone_pattern = Chem.MolFromSmarts("[CH3X4,CH2X4,CHX4]C(=O)[CH3X4,CH2X4,CHX4]")  # Ketone specific

    # Finding non-overlapping substructure matches in the molecule
    ketone_matches = mol.GetSubstructMatches(ketone_pattern, uniquify=True)

    # Extract unique center carbon atoms (carbonyl carbon atom) to avoid overlapping biases
    # Ensure extraction accurately represents distinct ketone occurrences
    unique_ketone_matches = {match[1] for match in ketone_matches}

    # Count unique ketone functionalities by core atom presence
    num_ketones = len(unique_ketone_matches)

    if num_ketones == 2:
        return True, "Molecule contains exactly two ketone functionalities"
    else:
        return False, f"Found {num_ketones} ketone functionalities, expected exactly 2"
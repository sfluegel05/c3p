"""
Classifies: CHEBI:27096 trichlorobenzene
"""
from rdkit import Chem

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule contains a trichlorobenzene core based on its SMILES string.
    The core is any benzene ring carrying three chloro substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a trichlorobenzene core, False otherwise.
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for any benzene ring with three Cl atoms
    trichloro_pattern = Chem.MolFromSmarts("c1cc(Cl)c(Cl)c(Cl)c1")
    
    # Check if the molecule contains a substructure that matches the pattern
    if mol.HasSubstructMatch(trichloro_pattern):
        return True, "Contains a benzene ring with exactly three chloro substituents"
    
    return False, "Does not contain a benzene ring with exactly three chloro substituents"
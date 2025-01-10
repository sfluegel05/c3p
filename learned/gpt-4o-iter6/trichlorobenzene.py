"""
Classifies: CHEBI:27096 trichlorobenzene
"""
from rdkit import Chem

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule contains a trichlorobenzene core based on its SMILES string.
    The core is any benzene ring carrying three chloro substituents at unspecified positions.

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

    # Define a SMARTS pattern for a benzene ring
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    
    # Check if the molecule contains a benzene ring
    if not mol.HasSubstructMatch(benzene_pattern):
        return False, "No benzene ring found"
    
    # Count chlorine atoms in the molecule
    chlorine_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "Cl")
    
    # Check for exactly three chlorine atoms in the molecule
    if chlorine_count != 3:
        return False, f"Contains {chlorine_count} Cl atoms, requires exactly 3"

    return True, "Contains a benzene ring with exactly three chloro substituents at unspecified positions"
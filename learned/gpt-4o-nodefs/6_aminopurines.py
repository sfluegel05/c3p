"""
Classifies: CHEBI:20706 6-aminopurines
"""
from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule belongs to the 6-aminopurines class based on its SMILES string.
    6-aminopurines are characterized by a purine ring with an amino group at the 6-position, possibly
    embedded in larger CoA-like structures.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str) - True if molecule is a 6-aminopurine, False otherwise with reason
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define an extended 6-aminopurine pattern
    # Considers an amino group attached to a purine-like ring possibly within a larger structure
    aminopurine_pattern = Chem.MolFromSmarts("Nc1ncnc2ncnc12")  # Amino purine structure

    # Check for 6-aminopurine substructure
    if mol.HasSubstructMatch(aminopurine_pattern):
        # Additionally check for larger contexts, like CoA associated structures.
        coa_like_pattern = Chem.MolFromSmarts("C(C)(COP(O)(=O)OP(O)(=O)")
        if mol.HasSubstructMatch(coa_like_pattern):
            return True, "Contains 6-aminopurine and CoA-like structures"

        return True, "Contains 6-aminopurine substructure"

    return False, "Does not contain 6-aminopurine substructure"
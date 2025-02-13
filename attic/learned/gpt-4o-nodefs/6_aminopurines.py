"""
Classifies: CHEBI:20706 6-aminopurines
"""
from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule belongs to the 6-aminopurines class based on its SMILES string.
    6-aminopurines are characterized by a purine ring with an amino group at the 6-position,
    potentially found within larger complex molecules like those associated with CoA.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str) - True if molecule is a 6-aminopurine, False otherwise with reason
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the 6-aminopurine structure itself
    # The pattern must match a purine structure with an amino group attached at the 6-position
    aminopurine_pattern = Chem.MolFromSmarts("Nc1ncnc2ncnc12")  # Amino group at position 6 of purine
    
    # Verify if 6-aminopurine substructure is present
    if not mol.HasSubstructMatch(aminopurine_pattern):
        return False, "Does not contain 6-aminopurine substructure"

    # Additional check for larger contexts or related structures
    # Ensure it doesn't mistakenly identify something else as the complex structure
    coa_like_pattern = Chem.MolFromSmarts("C[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n[cR1]c[nR1]c2n[cR1]nc2")
    
    # Check for complex structures associated commonly with 6-aminopurines
    if mol.HasSubstructMatch(coa_like_pattern):
        return True, "Contains 6-aminopurine and associated CoA-like structure"

    return True, "Contains 6-aminopurine substructure"
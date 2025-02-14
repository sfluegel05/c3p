"""
Classifies: CHEBI:27300 vitamin D
"""
from rdkit import Chem

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D compound based on its SMILES string.
    Vitamin D compounds are identified by their hydroxy seco-steroid backbone
    and specific structural features that are characteristic of this family of molecules.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a vitamin D compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Characteristic seco-steroid pattern observed in vitamin D
    seco_steroid_pattern = Chem.MolFromSmarts("C1CCC2=C(C1)C(C)(C)CCC3=CC=C4CCC2C34")
    if not mol.HasSubstructMatch(seco_steroid_pattern):
        return False, "No characteristic seco-steroid structure found"

    # Check for necessary hydroxyl groups, allowing for some variability
    hydroxyl_patterns = [
        Chem.MolFromSmarts("[C@H](O)"),  # Common chiral center with hydroxyl
        Chem.MolFromSmarts("[C@](C)(O)") # Another common hydroxyl pattern
    ]
    hydroxyl_count = 0
    for pattern in hydroxyl_patterns:
        if mol.HasSubstructMatch(pattern):
            hydroxyl_count += 1
    
    if hydroxyl_count < 1:
        return False, f"Insufficient hydroxyl group matches, found {hydroxyl_count}"

    # Identify long alkyl side chain typically part of vitamin D
    alkyl_chain_pattern = Chem.MolFromSmarts("C(C)CCC(C)CCC(O)C")
    if not mol.HasSubstructMatch(alkyl_chain_pattern):
        return False, "Missing typical alkyl side chain"

    return True, "Molecule matches key structural features of vitamin D"
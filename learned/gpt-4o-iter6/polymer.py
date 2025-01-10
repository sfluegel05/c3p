"""
Classifies: CHEBI:60027 polymer
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_polymer(smiles: str):
    """
    Determines if a molecule is a polymer based on its SMILES string.
    A polymer is characterized by repetitive structural units, often with high molecular weight.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polymer, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular weight - polymers typically have a high molecular weight
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for a typical polymer"

    # Look for repetitive patterns
    # Example: Repeating alkene units which are common in synthetic polymers
    repeating_unit_pattern = Chem.MolFromSmarts("C=C")
    repeating_matches = mol.GetSubstructMatches(repeating_unit_pattern)
    if len(repeating_matches) >= 3:
        return True, f"Contains {len(repeating_matches)} repeating units of an alkene"

    # Check for long carbon chains as an example of a polymer backbone
    carbon_chain_pattern = Chem.MolFromSmarts("C~C~C~C~C~C")  # Arbitrary chain length example
    if mol.HasSubstructMatch(carbon_chain_pattern):
        return True, "Contains long carbon chain, indicative of polymer backbone"

    # Optionally, count specific functional groups common in polymers
    # This is a simplified example; real polymer detection might require more sophisticated methods

    # If no indicative patterns are found, conclude it is not a polymer
    return False, "Does not match typical polymer characteristics"
"""
Classifies: CHEBI:60027 polymer
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_polymer(smiles: str):
    """
    Determines if a molecule is a polymer based on its SMILES string.
    A polymer is characterized by repetitive structural units, potentially high molecular weight, and 
    can have specific types of linkages.

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

    # Consider excluding entities with very low molecular weight first
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 200:  # Adjusted threshold to capture more polymers
        return False, "Molecular weight too low for a typical polymer"
    
    # Look for various repeating unit patterns
    repeating_unit_patterns = [
        Chem.MolFromSmarts("C=C"),  # Common repeating unit in vinyl-based polymers
        Chem.MolFromSmarts("C-C-C-C-C"),  # Example of a simple long carbon chain
        Chem.MolFromSmarts("O-C(=O)C"),  # Ester linkage, common in polyesters
        Chem.MolFromSmarts("C-O-C")  # Ether linkages, such as in polyethylene glycol
    ]
    
    found_repeating_unit = False
    for pattern in repeating_unit_patterns:
        repeating_matches = mol.GetSubstructMatches(pattern)
        if len(repeating_matches) >= 3:
            found_repeating_unit = True
            break

    if found_repeating_unit:
        return True, "Contains repeating units indicative of polymer structure"

    return False, "Does not match typical polymer characteristics"
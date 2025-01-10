"""
Classifies: CHEBI:67194 cannabinoid
"""
from rdkit import Chem

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    Cannabinoids generally feature a phenolic structure connected to 
    a potentially polyunsaturated hydrophobic side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cannabinoid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Define a basic cannabinoid pattern (example pattern covering phenolic structure and tail)
    # Note: This pattern is overly simplistic and may not correctly cover all cannabinoids
    phenolic_pattern = Chem.MolFromSmarts("c1cc(O)ccc1") # Phenolic ring
    if not mol.HasSubstructMatch(phenolic_pattern):
        return False, "No phenolic structure found, generally common in cannabinoids"

    # Check for a long hydrophobic chain often present in cannabinoids (possibly unsaturated)
    # Here we simply look for a series of CC repeating units suggesting a long chain
    hydrophobic_chain_pattern = Chem.MolFromSmarts("C(C)C")
    chain_matches = mol.GetSubstructMatches(hydrophobic_chain_pattern)
    if len(chain_matches) < 3:
        return False, "Insufficient hydrophobic chain, expected in many cannabinoids"
    
    # If both structural criteria are met, assume cannabinoid
    return True, "Contains a phenolic structure and a hydrophobic chain, typical of many cannabinoids"

# Note: This function is intended as a starting point. Cannabinoid classification often requires
# additional structural rules and data curation.
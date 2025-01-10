"""
Classifies: CHEBI:61384 sulfolipid
"""
from rdkit import Chem

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    A sulfolipid contains a sulfonic acid residue joined by a carbon-sulfur bond to a lipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Updated pattern for sulfonic acid to improve capture quality
    # Allow for varied attachment configurations (laughteral examples showed 
    # R-S(=O)(=O)[O-] where R indicates variations of chains)
    sulfonic_acid_pattern = Chem.MolFromSmarts("[S](=[O])(=[O])[O][C]")

    # Improved match ensuring split recognition of sulfonic moiety and assumed flexible R chain.
    if not mol.HasSubstructMatch(sulfonic_acid_pattern):
        return False, "Does not match a pattern of sulfonic acid properly joined to a carbon group"

    # Additional heuristic check for a broader lipid characteristic.
    # Check the overall presence of long carbon chains (implies a more general lipophilic character)
    carbon_chain_pattern = Chem.MolFromSmarts("C(-C)(-C)(-C)(-C)(-C)(-C)(-C)(-C)(-C)")

    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "Insufficient presence of a long carbon chain typically characteristic of lipids"
    
    return True, "Contains a sulfonic acid group connected by a carbon-sulfur bond to a lipid chain"
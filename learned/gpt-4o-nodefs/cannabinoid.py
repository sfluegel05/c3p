"""
Classifies: CHEBI:67194 cannabinoid
"""
from rdkit import Chem

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    Cannabinoids often feature complex ring structures, long unsaturated 
    carbon chains, and phenolic groups.

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
    
    # Define improved cannabinoid patterns
    # Pattern A: Phenolic ring typical in THC and CBD
    phenolic_pattern = Chem.MolFromSmarts("c1cc(O)ccc1")
    
    # Pattern B: Specific cannabinoid-like carbon tail length pattern (around 12-20 carbons)
    carbon_tail_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCC")  # Checking C12-C20
    
    # Pattern C: Basic tricyclic structure typical in few cannabinoids (classical cannabinol)
    tricyclic_structure_pattern = Chem.MolFromSmarts("C1=CC=C2C(=C1)C3=CC=CC=C3O2")  # Example pattern
    
    # Check any substructure matches
    phenolic_match = mol.HasSubstructMatch(phenolic_pattern)
    carbon_tail_match = mol.HasSubstructMatch(carbon_tail_pattern)
    tricyclic_match = mol.HasSubstructMatch(tricyclic_structure_pattern)
    
    # Implementing logic: Recognizing possible structural diversity in cannabinoids
    if phenolic_match or carbon_tail_match or tricyclic_match:
        reason = "Potential cannabinoid structure found with:"
        if phenolic_match:
            reason += " phenolic ring"
        if carbon_tail_match:
            if phenolic_match:
                reason += " and"
            reason += " long carbon chain"
        if tricyclic_match:
            if phenolic_match or carbon_tail_match:
                reason += " and"
            reason += " tricyclic structure"

        return True, reason

    return False, "No specific cannabinoid structure detected"
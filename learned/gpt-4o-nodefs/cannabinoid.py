"""
Classifies: CHEBI:67194 cannabinoid
"""
from rdkit import Chem

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    Cannabinoids often feature distinctive structural motifs such as a repeating alkyl 
    side chain with possible saturation, and sometimes a phenolic ring.

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
    
    # Define potential cannabinoid patterns
    # Pattern A: Phenolic ring common in classical cannabinoids (CBD, THC)
    phenolic_pattern = Chem.MolFromSmarts("c1cc(O)ccc1")
    
    # Pattern B: Long carbon chain often observed in cannabinoids (e.g., arachidonyl tail)
    hydrophobic_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCCC")  # Verify presence of a long C chain

    # Check any substructure matches
    phenolic_match = mol.HasSubstructMatch(phenolic_pattern)
    hydrophobic_match = mol.HasSubstructMatch(hydrophobic_chain_pattern)
    
    # Implementing logic: Recognize cannabinoids without phenolic groups but possessing alkyl chains
    if phenolic_match or hydrophobic_match:
        reason = "Potential cannabinoid structure found with:"
        if phenolic_match:
            reason += " phenolic ring"
        if hydrophobic_match:
            if phenolic_match:
                reason += " and"
            reason += " long carbon chain"

        return True, reason

    return False, "No specific cannabinoid structure detected"

# Note: This function takes a broader approach acknowledging possible cannabinoid diversity,
# including common structural motifs, but may still require refinement with further data.
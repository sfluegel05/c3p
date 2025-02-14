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
    
    # Match pattern for sulfonic acid attached to carbon indicating a sulfolipid
    sulfonic_lipid_pattern = Chem.MolFromSmarts("[C]-[S](=[O])(=[O])[O]")
    
    if not mol.HasSubstructMatch(sulfonic_lipid_pattern):
        return False, "No sulfonic acid group with a suitable carbon linkage to suggest a lipid structure found"

    # Check for long hydrocarbon chain for lipid property
    long_hydrocarbon_chain_pattern = Chem.MolFromSmarts("[C]-[S]([O])(=[O])=[O]C[CCCCCCCCCC]")  # Ensure presence of a longer chain
    if not mol.HasSubstructMatch(long_hydrocarbon_chain_pattern):
        return False, "No long carbon chain linked through a C-S bond suitable for lipid identification found"

    return True, "Contains a sulfonic acid group connected by a carbon-sulfur bond to a lipid chain"
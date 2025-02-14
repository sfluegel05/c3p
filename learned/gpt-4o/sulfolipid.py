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

    # Check for sulfonic acid group with connected carbon
    sulfonic_acid_lipid_pattern = Chem.MolFromSmarts("C-S(=O)(=O)-O")
    if not mol.HasSubstructMatch(sulfonic_acid_lipid_pattern):
        return False, "No sulfonic acid group with connected carbon (potential lipid linkage) found"

    # Verify carbon chain indicative of lipid presence
    carbon_chain_pattern = Chem.MolFromSmarts("C-[CX4,CX3]-[CX4,CX3,CX2,CX1]~[SX]-[OX2](=O)=O")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No appropriate carbon chain connected through C-S bond to a sulfonic group suggesting lipid"

    return True, "Contains a sulfonic acid group connected by a carbon-sulfur bond to a lipid"
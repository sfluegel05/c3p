"""
Classifies: CHEBI:63534 monoamine
"""
from rdkit import Chem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine is defined as an arylamino compound which contains one amino group connected
    to an aromatic ring by a two-carbon chain.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoamine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Updated SMARTS patterns for an amino group connected to an aromatic system
    # Primary amine SMARTS: N([H])(H)
    two_carbon_amino_pattern_primary = Chem.MolFromSmarts("[NX3;H2]-[CH2]-[CH2]-[c]") 
    # Include patterns for secondary amines
    two_carbon_amino_pattern_secondary = Chem.MolFromSmarts("[NX3;H1]-[CH2]-[CH2]-[c]") 

    # Check for two-carbon chain with amino group connected to aromatic system
    if mol.HasSubstructMatch(two_carbon_amino_pattern_primary):
        return True, "Contains primary amino group connected to aromatic ring by a two-carbon chain"
    elif mol.HasSubstructMatch(two_carbon_amino_pattern_secondary):
        return True, "Contains secondary amino group connected to aromatic ring by a two-carbon chain"

    return False, "Missing two-carbon chain with amino group connected to an aromatic ring"
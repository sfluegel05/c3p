"""
Classifies: CHEBI:63534 monoamine
"""
from rdkit import Chem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine typically has an amino group linked to an aromatic system
    by a two-carbon chain.

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

    # Look for a general aromatic pattern
    aromatic_pattern = Chem.MolFromSmarts("[a]") # any aromatic atom
    if not mol.HasSubstructMatch(aromatic_pattern):
        return False, "No aromatic system found"
    
    # Look for a generalized amine pattern linked through a two-carbon chain
    # [NX3;H2,H1,H0] matches any nitrogen (primary, secondary, tertiary amines)
    # connected via a two-carbon chain to the aromatic system
    amine_chain_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0]CC[a]")
    if not mol.HasSubstructMatch(amine_chain_pattern):
        return False, "No amine linked by a two-carbon chain to an aromatic system found"

    return True, "Contains an aromatic system with an amino group linked via a two-carbon chain"
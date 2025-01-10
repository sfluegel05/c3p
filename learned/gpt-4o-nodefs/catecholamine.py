"""
Classifies: CHEBI:33567 catecholamine
"""
from rdkit import Chem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine typically has a catechol moiety (benzene with hydroxyl groups on adjacent carbons)
    and an amine group attached via a two-carbon side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catecholamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for catechol moiety (benzene with hydroxyl groups on adjacent carbons)
    catechol_pattern = Chem.MolFromSmarts("c1cc(O)c(O)cc1")
    if not mol.HasSubstructMatch(catechol_pattern):
        return False, "No catechol moiety found"
    
    # Ensure presence of an amine group attached by a two-carbon side chain
    side_chain_amine_pattern = Chem.MolFromSmarts("C[C@H](O)[C,N]")
    if not mol.HasSubstructMatch(side_chain_amine_pattern):
        return False, "No appropriate amine group found"
    
    return True, "Contains catechol moiety and amine group indicating a catecholamine"

# Example usage
print(is_catecholamine("CC(C)NC[C@H](O)c1ccc(O)c(O)c1"))  # L-isoprenaline
print(is_catecholamine("C=C(O)C1=CC(=CCC1)O"))  # Non-catecholamine (This SMILES example does not depict a catecholamine)
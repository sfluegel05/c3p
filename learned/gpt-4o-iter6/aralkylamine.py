"""
Classifies: CHEBI:18000 aralkylamine
"""
from rdkit import Chem

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    An aralkylamine is defined as an alkylamine in which the alkyl group
    is substituted by an aromatic group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aralkylamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for a primary, secondary, tertiary, or quaternary amine, excluding amides
    amine_pattern = Chem.MolFromSmarts("[NX3+0,NX4+;H2,H1,H0;!$([NX3][C]=O)]")
    if not mol.HasSubstructMatch(amine_pattern):
        return False, "No amine group found"
    
    # Updated pattern: Aromatic group must be connected through non-ring atoms to an amine group
    # This accounts for aromatic carbon attachment to alkyl, then to amine
    aralkylamine_pattern = Chem.MolFromSmarts("[#6][c]~[*][CH2][NX3+0,NX4+;H2,H1,H0;!$([NX3+0][C]=O)]")
    if not mol.HasSubstructMatch(aralkylamine_pattern):
        return False, "No aralkylamine feature found"

    return True, "Contains alkylamine group with aromatic substitution"
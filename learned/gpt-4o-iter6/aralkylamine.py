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
    
    # Look for a primary, secondary, or tertiary amine not part of an amide
    amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$([NX3][C]=O)]")
    if not mol.HasSubstructMatch(amine_pattern):
        return False, "No amine group found"
    
    # Check for an aromatic group connected to an alkyl chain (with an amine)
    aralkylamine_pattern = Chem.MolFromSmarts("[c]~[CX4,CH,CH2]~[NX3;H2,H1,H0;!$([NX3][C]=O)]") 
    if not mol.HasSubstructMatch(aralkylamine_pattern):
        return False, "No aralkylamine feature found"

    return True, "Contains alkylamine group with aromatic substitution"
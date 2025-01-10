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
    
    # Look for a nitrogen with at least one alkyl group not part of an amide
    amine_pattern = Chem.MolFromSmarts("[NX3;!$(NC=O)]")
    amine_matches = mol.GetSubstructMatches(amine_pattern)
    if not amine_matches:
        return False, "No amine group found"
    
    # Look for alkyl attached to aromatic carbon adjacent to amine
    aralkyl_pattern = Chem.MolFromSmarts("[NX3;!$(NC=O)]~[CX4]~a")  # Amine connected to an alkyl which is connected to an aromatic
    if not mol.HasSubstructMatch(aralkyl_pattern):
        return False, "Alkylamine without aromatic substitution found"
    
    return True, "Contains alkylamine group with aromatic substitution"
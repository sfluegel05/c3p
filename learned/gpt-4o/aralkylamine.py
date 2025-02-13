"""
Classifies: CHEBI:18000 aralkylamine
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    An aralkylamine is defined as an alkylamine in which the alkyl group is substituted by an aromatic group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aralkylamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define patterns for aromatic and amine groups
    aromatic_pattern = Chem.MolFromSmarts("a")
    amine_patterns = [
        Chem.MolFromSmarts("[NX3;H2]"),  # Primary amine
        Chem.MolFromSmarts("[NX3;H1]"),  # Secondary amine
        Chem.MolFromSmarts("[NX3;H0]"),  # Tertiary amine
    ]
    
    # Check for at least one aromatic ring
    if not mol.HasSubstructMatch(aromatic_pattern):
        return False, "No aromatic group found"
    
    # Check for presence of amine groups
    for amine_pattern in amine_patterns:
        if mol.HasSubstructMatch(amine_pattern):
            # Now verify the linkage of the aromatic group via alkyl chain
            # We'll look for a generic pattern: Aromatic - Alkyl - NHx
            aralkylamine_pattern = Chem.MolFromSmarts("a-[CX4]-[NX3]")
            if mol.HasSubstructMatch(aralkylamine_pattern):
                return True, "Contains an alkylamine with aromatic substitution"
            else:
                return False, "Aromatic group not attached through an alkyl to amine"
    
    return False, "No amine group found"
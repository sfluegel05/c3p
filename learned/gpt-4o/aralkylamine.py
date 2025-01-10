"""
Classifies: CHEBI:18000 aralkylamine
"""
from rdkit import Chem

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    An aralkylamine is defined as an alkylamine in which the alkyl group is substituted by an aromatic group.

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
    
    # Look for primary or secondary amine pattern (N not part of carbonyl)
    amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")
    if not mol.HasSubstructMatch(amine_pattern):
        return False, "No amine group found"
    
    # Look for aromatic group pattern
    aromatic_pattern = Chem.MolFromSmarts("a")
    if not mol.HasSubstructMatch(aromatic_pattern):
        return False, "No aromatic group found"

    # Look for the linkage pattern: Aromatic -> Aliphatic Carbon -> Amine
    linkage_pattern = Chem.MolFromSmarts("a-[CX4]-[NX3;H2,H1;!$(NC=O)]")
    if not mol.HasSubstructMatch(linkage_pattern):
        return False, "No linkage between aromatic group and amine via an alkyl chain found"

    return True, "Contains an amine group linked to an aromatic group by an alkyl chain"
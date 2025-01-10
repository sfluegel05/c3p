"""
Classifies: CHEBI:18000 aralkylamine
"""
from rdkit import Chem
from rdkit.Chem import rdpatterns

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
    
    # Look for amine pattern (nitrogen not in carbonyls and bonded to alkyl)
    amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")  # Primary, secondary amine
    if not mol.HasSubstructMatch(amine_pattern):
        return False, "No amine group found"
    
    # Look for aromatic ring
    aromatic_pattern = Chem.MolFromSmarts("a")  # Any aromatic atom
    if not mol.HasSubstructMatch(aromatic_pattern):
        return False, "No aromatic group found"
    
    # Verify connection through a non-aromatic alkyl chain
    alkyl_amino_pattern = Chem.MolFromSmarts("[CX4][NX3;H2,H1;!$(NC=O)]")
    # Look for linkage: Aromatic -> Aliphatic C -> N (Amine)
    combined_pattern = Chem.MolFromSmarts("*a-[CX4]-[NX3;H2,H1]") 
    
    if not mol.HasSubstructMatch(combined_pattern):
        return False, "No alkyl chain linking amine and aromatic group"

    return True, "Contains an amine group linked to an aromatic group by an alkyl chain"
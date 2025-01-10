"""
Classifies: CHEBI:35785 sphingoid
"""
from rdkit import Chem

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    Sphingoids consist of a long chained aliphatic component with hydroxyl and
    amino functionalities, along with possible unsaturations.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a sphingoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for a long aliphatic chain, allowing unsaturation and variations
    long_chain_pattern = Chem.MolFromSmarts("[C]([C])[C]([C])([C])[C]([C])([C])[C]([C])([C])[C]([C])([C])[C]")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "Does not have the correct patterns for a long aliphatic carbon chain"

    # Detecting primary or secondary amino group attached to carbon
    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1][CX4]")
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "Missing amino group"
    
    # Detect hydroxyl groups connected to carbon
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H][CX4]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Missing hydroxyl group"

    # Optional: Detect presence of double bonds in carbon chains
    double_bond_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    double_bonds_present = mol.HasSubstructMatch(double_bond_pattern)

    # Provide appropriate reasoning based on presence of unsaturation
    if double_bonds_present:
        return True, "Contains aliphatic chain with hydroxyl and amino groups, allowing for unsaturation; likely a sphingoid."
    else:
        return True, "Contains aliphatic chain with hydroxyl and amino groups even if fully saturated; likely a sphingoid."
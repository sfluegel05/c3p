"""
Classifies: CHEBI:35785 sphingoid
"""
from rdkit import Chem

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    Sphingoids consist of a long aliphatic chain, typically with hydroxyl and
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

    # Looser pattern for a long aliphatic chain, considering variations
    long_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCC")  # Example: matches at least 10 carbon atoms in any configuration
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "Does not have the correct patterns for a long aliphatic carbon chain"

    # Detecting primary, secondary, or tertiary amino group attached to carbon
    amino_pattern = Chem.MolFromSmarts("[NX3][CX4;!$(C=O)]")
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "Missing amino group"
    
    # Detect hydroxyl groups connected to carbon
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2][CX4]")  # More general pattern for detection
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Missing hydroxyl group"

    # Presence of double bonds in carbon chains
    double_bond_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    double_bonds_present = mol.HasSubstructMatch(double_bond_pattern)

    # Provide appropriate reasoning based on presence of unsaturation
    if double_bonds_present:
        return True, "Contains long aliphatic chain with hydroxyl and amino groups, possibly unsaturated; likely a sphingoid."
    else:
        return True, "Contains long aliphatic chain with hydroxyl and amino groups, likely a sphingoid."

# Note: The efficacy of the detection may require adjustment of patterns based 
# on specific sphingoid variations not fully captured by general patterns.
"""
Classifies: CHEBI:35785 sphingoid
"""
from rdkit import Chem

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    Sphingoids are defined primarily by the presence of a long aliphatic chain with
    one or more hydroxyl groups and an amino group, along with possible unsaturations.

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
    
    # Check for long aliphatic chain
    carbon_chain_pattern = Chem.MolFromSmarts("[CCCCCCCCC;R0]")  # Approx >6 carbons in a row
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "Does not have a long aliphatic carbon chain"

    # Check for hydroxyl group(s)
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4](O)")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 1:
        return False, "Missing hydroxyl group"

    # Check for an amino group
    amino_pattern = Chem.MolFromSmarts("[CX4,NX3H2]")
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    if len(amino_matches) < 1:
        return False, "Missing amino group"
    
    # Check for double bonds (optional for unsaturation)
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "Lacks unsaturation, might still be a sphingoid"

    # Verify proper positioning of functional groups relative to chain length
    c_count = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if c_count < 14:
        return False, "Aliphatic chain too short for typical sphingoid"
    
    return True, "Contains long aliphatic chain with hydroxyl and amino groups; could be a sphingoid"
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

    # Updated pattern for detecting aliphatic chains, considering both saturated and unsaturated
    long_chain_pattern = Chem.MolFromSmarts("[CH2]1[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]1")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "Does not have the correct patterns for a long aliphatic carbon chain"

    # Broadening criteria for amino groups, considering position and type
    primitive_amino_patterns = ["[CX3][NX3;H2,H1;!$(NC=O)]", "[NH2,NH1]"]
    amino_present = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in primitive_amino_patterns)
    if not amino_present:
        return False, "Missing amino group"
    
    # Broaden hydroxyl group detection - check for adjacent to chiral centers
    hydroxyl_patterns = ["[CX3]([OH])", "[O][CX3;H2]"]
    hydroxyl_present = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in hydroxyl_patterns)
    if not hydroxyl_present:
        return False, "Missing hydroxyl group"

    # Sphingoids can optionally include double bonds (usually one or more)
    double_bond_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    if mol.HasSubstructMatch(double_bond_pattern):
        return True, "Contains aliphatic chain with hydroxyl and amino groups, allowing unsaturation; could be a sphingoid."
    else:
        return True, "Contains aliphatic chain with hydroxyl and amino groups even if fully saturated; it could be a sphingoid."
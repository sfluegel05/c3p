"""
Classifies: CHEBI:17636 sphingomyelin d18:1
"""
from rdkit import Chem

def is_sphingomyelin_d18_1(smiles: str):
    """
    Determines if a molecule is a sphingomyelin d18:1 based on its SMILES string.
    'Sphingomyelin d18:1' is defined as having sphingosine as the sphingoid component.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingomyelin d18:1, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More general sphingosine backbone pattern
    # Allow variations in stereochemistry and chain length near the double bond
    sphingosine_pattern = Chem.MolFromSmarts("O[C@H](C)[C@H](N)C=C")  # Simplified sphingosine core
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone (d18:1) pattern found"

    # Check for the phosphocholine group
    phosphocholine_pattern = Chem.MolFromSmarts("COP([O-])([O-])=O")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine moiety found"
    
    # Look for at least one amide bond (usually, an acyl group is attached to the nitrogen)
    amide_pattern = Chem.MolFromSmarts("NC(=O)")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage found"

    return True, "Molecule has sphingosine backbone and phosphocholine group"
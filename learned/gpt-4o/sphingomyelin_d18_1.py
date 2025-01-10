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

    # Advanced sphingosine backbone pattern with specific stereochemistry
    sphingosine_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](CC=C)NC1")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone (d18:1) pattern found"

    # Check for the phosphocholine group
    phosphocholine_pattern = Chem.MolFromSmarts("O=P(OCC[N+](C)(C)C)(O)O")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine moiety found"
    
    # Look for at least one amide bond (usually, an acyl group is attached to the nitrogen)
    amide_pattern = Chem.MolFromSmarts("NC(=O)")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage found"

    return True, "Molecule has sphingosine backbone, phosphocholine group, and amide linkage"
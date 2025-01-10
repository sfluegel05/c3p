"""
Classifies: CHEBI:17636 sphingomyelin d18:1
"""
from rdkit import Chem

def is_sphingomyelin_d18_1(smiles: str):
    """
    Determines if a molecule is a sphingomyelin d18:1 based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a sphingomyelin d18:1, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Correct sphingosine backbone pattern
    # Pattern: C1[C@H](O)[C@H](N)[C]=[C]1
    sphingosine_pattern = Chem.MolFromSmarts("C[C@@H](O)C[N]C=C")
    if sphingosine_pattern is None or not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone (sphingoid base) found"
    
    # Look for amide bond -C(=O)N, linking fatty acid to sphingosine
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if amide_pattern is None or not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group (linking sphingosine and fatty acid) found"
    
    # Ensure phosphocholine headgroup presence -COP(=O)([O-])OCC[N+](C)(C)C
    phosphocholine_pattern = Chem.MolFromSmarts("COP(=O)([O-])OCC[N+](C)(C)C")
    if phosphocholine_pattern is None or not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine headgroup found"

    return True, "Contains sphingosine backbone with fatty acid and phosphocholine headgroup"
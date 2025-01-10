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

    # Look for sphingosine backbone (sphingoid base)
    # Pattern: C1[C@@H](O)[C@H](COP(=O)([O-]))NC(=O)
    sphingosine_pattern = Chem.MolFromSmarts("[C@@H](O)[C@H](N)[C=C][C]O")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone (sphingoid base) found"
    
    # Look for amide bond -C(=O)N
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group (linking sphingosine and fatty acid) found"
    
    # Look for phosphocholine headgroup -COP(=O)([O-])OCC[N+](C)(C)C
    phosphocholine_pattern = Chem.MolFromSmarts("COP(=O)([O-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine headgroup found"

    return True, "Contains sphingosine backbone with fatty acid and phosphocholine headgroup"
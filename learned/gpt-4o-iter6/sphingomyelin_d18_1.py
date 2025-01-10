"""
Classifies: CHEBI:17636 sphingomyelin d18:1
"""
from rdkit import Chem

def is_sphingomyelin_d18_1(smiles: str):
    """
    Determines if a molecule is a sphingomyelin d18:1 based on its SMILES string.
    
    A sphingomyelin d18:1 contains a sphingosine backbone, an amide-linked fatty acid,
    and a phosphocholine headgroup.
    
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

    # Sphingosine backbone (d18:1) pattern with the right stereochemistry
    sphingosine_pattern = Chem.MolFromSmarts("[C@@H](O)[C@@H](NC(=O))CCC/C=C\\")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone (sphingoid base) found matching d18:1 structure"
    
    # Amide bond -C(=O)N pattern
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group (linking sphingosine and fatty acid) found"
    
    # Phosphocholine headgroup pattern
    phosphocholine_pattern = Chem.MolFromSmarts("COP(=O)([O-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine headgroup found"

    return True, "Contains sphingosine backbone (d18:1) with fatty acid and phosphocholine headgroup"
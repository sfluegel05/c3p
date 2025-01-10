"""
Classifies: CHEBI:17636 sphingomyelin d18:1
"""
from rdkit import Chem

def is_sphingomyelin_d18_1(smiles: str):
    """
    Determines if a molecule is a sphingomyelin d18:1 based on its SMILES string.
    
    A sphingomyelin d18:1 contains a sphingosine backbone, which is characterized by 18 carbon atoms,
    one double bond, specific chirality, an amide-linked fatty acid, and a phosphocholine headgroup.
    
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

    # Define refined sphingosine backbone (d18:1) with explicit stereochemistry and double bond
    sphingosine_pattern = Chem.MolFromSmarts("[C@H](O)C[C@@H](NC(=O))C=CC")
    if sphingosine_pattern is None:
        return None, None
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone (sphingoid base) found with correct d18:1 structure"
    
    # Check for length of the chain
    chain_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCCC")  # 15 or more carbons in a chain
    if chain_pattern is None:
        return None, None
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No long hydrocarbon chain matching sphingosine characteristic found"
    
    # Amide bond -C(=O)N pattern
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if amide_pattern is None:
        return None, None
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group (linking sphingosine and fatty acid) found"
    
    # Phosphocholine headgroup pattern
    phosphocholine_pattern = Chem.MolFromSmarts("COP(=O)([O-])OCC[N+](C)(C)C")
    if phosphocholine_pattern is None:
        return None, None
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine headgroup found"

    return True, "Contains sphingosine backbone (d18:1) with fatty acid and phosphocholine headgroup"
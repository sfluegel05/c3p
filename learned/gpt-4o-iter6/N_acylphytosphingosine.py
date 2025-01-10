"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
from rdkit import Chem

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    An N-acylphytosphingosine is a ceramide with a phytosphingosine backbone 
    having a fatty acyl group attached to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylphytosphingosine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Phytosphingosine backbone pattern (simplified for structural recognition)
    phytosphingosine_pattern = Chem.MolFromSmarts("C[C@@H](O)[C@H](O)N")
    if not mol.HasSubstructMatch(phytosphingosine_pattern):
        return False, "No phytosphingosine backbone found"

    # Fatty acyl amide linkage pattern - R-C(=O)-N
    acyl_amide_pattern = Chem.MolFromSmarts("C(=O)[NH1]")
    if not mol.HasSubstructMatch(acyl_amide_pattern):
        return False, "No fatty acyl amide group found attached to nitrogen"
    
    # Additional considerations could involve ensuring long carbon chain lengths
    # but are omitted for simplicity in this implementation.

    return True, "Molecule contains a phytosphingosine backbone with a fatty acyl amide group"
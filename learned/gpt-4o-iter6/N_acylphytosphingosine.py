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
    
    # Improved pattern for phytosphingosine backbone
    phytosphingosine_pattern = Chem.MolFromSmarts("[C@H](CO)NC(=O)[C@H](O)[C@H](O)C")
    if not mol.HasSubstructMatch(phytosphingosine_pattern):
        return False, "No phytosphingosine backbone found"
    
    # Check for presence of an amide-bound long chain fatty acyl group
    # We will consider a simple carbon chain connected through an amide linkage
    acyl_amide_pattern = Chem.MolFromSmarts("C(=O)N[C@H]")
    if not mol.HasSubstructMatch(acyl_amide_pattern):
        return False, "No fatty acyl group bound via amide linkage"

    # Considering the biochemical diversity, ensuring acyl chain has sufficient carbon length typical of lipids
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 20:
        return False, "The molecule does not have enough carbon atoms for an N-acylphytosphingosine"
    
    return True, "The structure matches an N-acylphytosphingosine, with a phytosphingosine backbone and fatty acyl group"
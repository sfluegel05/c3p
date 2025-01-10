"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    A N-acylphytosphingosine has a phytosphingosine backbone with a fatty acyl group attached to nitrogen.

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

    # Look for phytosphingosine backbone: C-C-C-C with multiple hydroxyls and an amino group
    # Note: Phytosphingosine pattern with flexibility in chain length and hydroxyl positions.
    phytosphingosine_pattern = Chem.MolFromSmarts("CCCC(CO)NC([C@@H](CO)O)O")
    if not mol.HasSubstructMatch(phytosphingosine_pattern):
        return False, "No phytosphingosine backbone found"
    
    # Look for fatty acyl group structure: long carbon chain ending in an amide linkage
    fatty_acyl_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3][CX4,CX3]")
    if not mol.HasSubstructMatch(fatty_acyl_pattern):
        return False, "No fatty acyl group attached to nitrogen found"

    return True, "Contains phytosphingosine backbone with fatty acyl group attached to nitrogen"
"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
from rdkit import Chem

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    An N-acylphytosphingosine has a phytosphingosine backbone with a fatty acyl group attached to nitrogen.

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

    # Refined pattern for a phytosphingosine backbone:
    phytosphingosine_pattern = Chem.MolFromSmarts("[#6]-[#6]-[#6]-[#6]-[#8]-[#6]-[#8]-[#6]-[#7]-[#6]-[#8]")
    if not mol.HasSubstructMatch(phytosphingosine_pattern):
        return False, "No phytosphingosine backbone found"

    # Look for the fatty acyl group: strict amide linkage to the nitrogen
    fatty_acyl_pattern = Chem.MolFromSmarts("[#6](=O)[#7]-[#6]")
    if not mol.HasSubstructMatch(fatty_acyl_pattern):
        return False, "No fatty acyl group attached to nitrogen found"

    # Ensure the presence of sufficient hydrophobic chain length
    carbon_atoms = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    if carbon_atoms < 30:  # Typical length for N-acylphytosphingosines including both parts
        return False, f"Overall carbon chain length {carbon_atoms} is too short for N-acylphytosphingosine"

    return True, "Contains a phytosphingosine backbone with a fatty acyl group attached to nitrogen"
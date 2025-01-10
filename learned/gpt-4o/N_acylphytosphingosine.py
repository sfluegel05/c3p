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

    # Improved pattern for a phytosphingosine backbone:
    # A long chain (typically C18) with an amino group and hydroxyls.
    # Let's assume a simple carbon N-OH backbone with chirality checks.
    phytosphingosine_pattern = Chem.MolFromSmarts("[#6]-[#6]-[#6]-[#6]-[#8]-[#6]-[#6]-[#6]-[#7]([#6])-[#6]-[#6]-[#8]-[#6]-[#8]")
    if not mol.HasSubstructMatch(phytosphingosine_pattern):
        return False, "No phytosphingosine backbone found"

    # Look for the fatty acyl group: amide linkage to the nitrogen (as seen in CN)
    fatty_acyl_pattern = Chem.MolFromSmarts("[#6](=O)[#7]-[#6]")
    if not mol.HasSubstructMatch(fatty_acyl_pattern):
        return False, "No fatty acyl group attached to nitrogen found"

    # Ensure the presence of long hydrophobic chains (indicative of phytosphingosine with N-acyl)
    c_chain_length = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    if c_chain_length < 20:
        return False, f"Carbon chain length {c_chain_length} is too short for N-acylphytosphingosine"

    return True, "Contains a phytosphingosine backbone with a fatty acyl group attached to nitrogen"
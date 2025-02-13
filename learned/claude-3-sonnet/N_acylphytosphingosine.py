"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
"""
Classifies: CHEBI:84670 N-acylphytosphingosine

A ceramide that is phytosphingosine having a fatty acyl group attached to the nitrogen.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.

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
    
    # Look for phytosphingosine backbone pattern with C-4 hydroxyl group
    phytosphingosine_pattern = Chem.MolFromSmarts("[NH1][CX4]([CH2][CH2][CH2][CH2][CH2])[CX4]([CH1]([OH1])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2])[CH2][OH1]")
    phytosphingosine_match = mol.GetSubstructMatches(phytosphingosine_pattern)
    if not phytosphingosine_match:
        return False, "Phytosphingosine backbone not found"
    
    # Look for acyl chain attached to N atom with optional C-2 hydroxyl group
    acyl_pattern = Chem.MolFromSmarts("[NH1][CX3](=[OX1])[CX3]([OH1])?[CX3]~[CX3]")
    acyl_match = mol.GetSubstructMatches(acyl_pattern)
    if not acyl_match:
        return False, "No acyl group attached to nitrogen"
    
    # Look for common substituents like galactose
    galactose_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)O)O)O)O")
    galactose_match = mol.GetSubstructMatches(galactose_pattern)
    
    return True, "Contains phytosphingosine backbone with acyl group attached to nitrogen" + (", and galactose substituent" if galactose_match else "")
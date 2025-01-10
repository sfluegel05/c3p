"""
Classifies: CHEBI:67197 endocannabinoid
"""
from rdkit import Chem

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an endocannabinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of long carbon chains with double bonds
    long_chain_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(long_chain_pattern)
    if len(double_bond_matches) < 2:
        return False, "Fewer than 2 double bonds found, not a characteristic polyunsaturated chain."

    # Check for an ethanolamine group (NCCO or NCC(=O))
    ethanolamine_pattern = Chem.MolFromSmarts("NCCO")
    amide_pattern = Chem.MolFromSmarts("NCC(=O)")
    if not (mol.HasSubstructMatch(ethanolamine_pattern) or mol.HasSubstructMatch(amide_pattern)):
        return False, "No ethanolamine or amide group found."

    # Check for ether or ester linkages for diverse endocannabinoid structures
    ether_pattern = Chem.MolFromSmarts("OCC")
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    if not (mol.HasSubstructMatch(ether_pattern) or mol.HasSubstructMatch(ester_pattern)):
        return False, "No ether or ester linkages found."

    return True, "Matches characteristics of known endocannabinoids."
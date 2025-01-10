"""
Classifies: CHEBI:29256 thiol
"""
from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol is an organosulfur compound in which a thiol group (-SH) or a thioether group (-S-R) is attached to a carbon atom of any aliphatic or aromatic moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Exclude peptides by checking for peptide bonds (amide bonds)
    peptide_bond_smarts = '[NX3][CX3](=O)'
    peptide_bond_pattern = Chem.MolFromSmarts(peptide_bond_smarts)
    if mol.HasSubstructMatch(peptide_bond_pattern):
        return False, "Molecule appears to be a peptide containing sulfur"

    # Exclude sulfur atoms bonded to oxygen or nitrogen (sulfoxides, sulfones, sulfates, sulfonamides)
    excluded_sulfur_smarts = '[#16;D2](=[O,N])[O,N]'
    excluded_sulfur_pattern = Chem.MolFromSmarts(excluded_sulfur_smarts)
    if mol.HasSubstructMatch(excluded_sulfur_pattern):
        return False, "Sulfur bonded to oxygen or nitrogen (possible sulfoxide, sulfone, sulfate, or sulfonamide)"

    # Define SMARTS pattern for sulfur attached to carbon via single bonds, not attached to oxygen or nitrogen
    organosulfur_smarts = '[#16;D2;!$(S(~O)~O);!$(S(~O)~N);!$(S(~N)~N)](-[#6])-[#6]'
    organosulfur_pattern = Chem.MolFromSmarts(organosulfur_smarts)
    if mol.HasSubstructMatch(organosulfur_pattern):
        return True, "Contains sulfur attached to carbon atom"

    # Also check for thiol groups (-SH attached to carbon)
    thiol_smarts = '[#16;H1]-[#6]'
    thiol_pattern = Chem.MolFromSmarts(thiol_smarts)
    if mol.HasSubstructMatch(thiol_pattern):
        return True, "Contains thiol group (-SH) attached to carbon atom"

    return False, "No sulfur attached to carbon atom found"
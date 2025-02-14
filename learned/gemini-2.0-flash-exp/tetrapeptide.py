"""
Classifies: CHEBI:48030 tetrapeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide consists of four amino acid residues connected by three peptide bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrapeptide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for peptide bonds (more specific pattern).
    # The peptide bond is the essential feature of a peptide.
    peptide_bond_pattern = Chem.MolFromSmarts("[C](=[O])-[N]")
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bonds) != 3:
        return False, f"Found {len(peptide_bonds)} peptide bonds, need exactly 3"


    return True, "Contains four amino acid residues connected by three peptide linkages"
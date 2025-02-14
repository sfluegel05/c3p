"""
Classifies: CHEBI:46761 dipeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide is a molecule consisting of two amino acid residues connected by a peptide bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dipeptide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a more specific peptide bond pattern.
    # This pattern looks for a C(=O)-N fragment between two alpha carbons.
    peptide_bond_pattern = Chem.MolFromSmarts("[CX4][CX3](=[OX1])N[CX4]")
    
    peptide_matches = mol.GetSubstructMatches(peptide_bond_pattern)


    if len(peptide_matches) != 1:
        return False, f"Molecule has {len(peptide_matches)} peptide bonds, a dipeptide should have exactly 1"


    return True, "Contains two amino acid residues connected by a single peptide bond"
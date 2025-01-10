"""
Classifies: CHEBI:48030 tetrapeptide
"""
from rdkit import Chem

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide must contain four amino-acid residues connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrapeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Search for amide linkages (peptide bonds) C(=O)N in the structure
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)

    # Count peptide (amide) linkages, should be 3 for a tetrapeptide
    number_of_amides = len(amide_matches)
    if number_of_amides == 3:
        return True, "Contains four amino-acid residues connected by peptide bonds"
    
    return False, f"Contains {number_of_amides+1} amino-acid residues; expected exactly 4"
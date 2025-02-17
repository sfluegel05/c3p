"""
Classifies: CHEBI:35358 sulfonamide
"""
"""
Classifies: Sulfonamide 
Definition: An amide of a sulfonic acid RS(=O)2NR'2, which features the S(=O)(=O)N substructure.
"""

from rdkit import Chem

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.
    A sulfonamide is defined as an amide of a sulfonic acid, having the substructure RS(=O)2NR'2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains a sulfonamide group, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Create a SMARTS pattern for the sulfonamide group.
    # The pattern S(=O)(=O)N corresponds to an S atom (sulfur) double bonded to two oxygens and bound to a nitrogen.
    sulfonamide_pattern = Chem.MolFromSmarts("S(=O)(=O)N")
    
    # Check whether the molecule contains the sulfonamide substructure.
    if mol.HasSubstructMatch(sulfonamide_pattern):
        return True, "Molecule contains the sulfonamide functional group S(=O)(=O)N"
    else:
        return False, "Molecule does not contain a sulfonamide group S(=O)(=O)N"
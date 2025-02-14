"""
Classifies: CHEBI:35358 sulfonamide
"""
from rdkit import Chem

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.
    A sulfonamide is an amide of a sulfonic acid RS(=O)2NR'2 where R can vary.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a sulfonamide, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a general sulfonamide
    # This pattern includes the core sulfonamide connectivity: S(=O)(=O)N
    # We allow for connectivity indicating R groups on both S and N
    sulfonamide_pattern = Chem.MolFromSmarts("S(=O)(=O)[N]")

    # Check for the presence of the sulfonamide functional group
    if mol.HasSubstructMatch(sulfonamide_pattern):
        return True, "Contains the sulfonamide functional group S(=O)(=O)N"
    else:
        return False, "Does not contain the sulfonamide functional group S(=O)(=O)N"

# Example usage
smiles_example = "COCCCNS(=O)(=O)C1=C(C(=C(C=C1)OC)Cl)Cl"
result, reason = is_sulfonamide(smiles_example)
print(result, reason)
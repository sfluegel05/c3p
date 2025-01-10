"""
Classifies: CHEBI:36916 cation
"""
from rdkit import Chem

def is_cation(smiles: str):
    """
    Determines if a molecule is a cation based on its SMILES string.
    A cation is defined as a monoatomic or polyatomic species
    having one or more elementary charges of the proton.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cation, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for positive charges
    positive_count = sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() > 0)
    
    if positive_count > 0:
        return True, f"Molecule has {positive_count} positive charge(s)"
    else:
        return False, "No positive charges found"
"""
Classifies: CHEBI:36916 cation
"""
"""
Classifies: Cations (species with one or more positive charges)
"""
from rdkit import Chem

def is_cation(smiles: str):
    """
    Determines if a molecule is a cation based on its SMILES string.
    A cation is a species with one or more positive charges.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cation, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check each atom for positive charge
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() > 0:
            return True, f"Atom {atom.GetIdx()+1} has a positive charge"
    
    # Check for overall molecular charge (uncommon in SMILES but possible)
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if total_charge > 0:
        return True, f"Molecular charge is +{total_charge}"
    
    return False, "No positive charges found"
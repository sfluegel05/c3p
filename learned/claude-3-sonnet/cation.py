"""
Classifies: CHEBI:36916 cation
"""
"""
Classifies: CHEBI:36916 cation
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_cation(smiles: str):
    """
    Determines if a molecule is a cation based on its SMILES string.
    A cation is a monoatomic or polyatomic species having one or more 
    elementary charges of the proton.

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

    # Calculate formal charges on each atom
    positive_charges = 0
    negative_charges = 0
    positive_atoms = []
    
    for atom in mol.GetAtoms():
        formal_charge = atom.GetFormalCharge()
        if formal_charge > 0:
            positive_charges += formal_charge
            positive_atoms.append(f"{atom.GetSymbol()}{atom.GetIdx()+1}")
        elif formal_charge < 0:
            negative_charges += abs(formal_charge)
    
    # Calculate net charge
    net_charge = positive_charges - negative_charges
    
    # If no charges present
    if positive_charges == 0:
        return False, "No positive charges found"
    
    # If net charge is not positive
    if net_charge <= 0:
        return False, f"Net charge ({net_charge}) is not positive"
        
    # Success case
    charge_locations = ", ".join(positive_atoms)
    return True, f"Found net positive charge ({net_charge}+) with positive charges on: {charge_locations}"
"""
Classifies: CHEBI:33521 metal atom
"""
from rdkit import Chem

def is_metal_atom(smiles: str):
    """
    Determines if a SMILES string represents a single metal atom.
    Validates that the structure is a single atom of a metallic element.

    Args:
        smiles (str): Input SMILES string

    Returns:
        bool: True if metal atom, False otherwise
        str: Reason for decision
    """
    mol = Chem.MolFromSmiles(smiles)
    
    # Check if SMILES is valid
    if not mol:
        return False, "Invalid SMILES"
    
    # Must have exactly one atom
    if mol.GetNumAtoms() != 1:
        return False, "Not a single atom"
    
    atom = mol.GetAtomWithIdx(0)
    
    # Check if the element is a metal
    if not atom.GetIsAromatic() and atom.GetAtomicNum() in Chem.GetPeriodicTable().GetMetals():
        # Check for charge - metals in elemental form should have 0 charge
        if atom.GetFormalCharge() == 0:
            return True, "Single metallic atom in elemental state"
        else:
            return False, "Charged metal species"
    
    return False, "Not a metallic element"
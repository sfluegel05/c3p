"""
Classifies: CHEBI:33521 metal atom
"""
from rdkit import Chem

def is_metal_atom(smiles: str):
    """
    Determines if a SMILES string represents a single metal atom in its elemental state.
    Validates that the structure is a single atom of a metallic element with 0 charge.

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
    
    # Check if the element is a metal using RDKit's built-in metal detection
    if not atom.GetIsMetal():
        return False, "Atom is not a metallic element"
    
    # Check for neutral charge (elemental form)
    if atom.GetFormalCharge() != 0:
        return False, "Charged metal species"
    
    # Verify no radical electrons (though rare in standard representations)
    if atom.GetNumRadicalElectrons() != 0:
        return False, "Radical present in metal atom"
    
    return True, "Single metallic atom in elemental state"
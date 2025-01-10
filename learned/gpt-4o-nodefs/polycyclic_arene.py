"""
Classifies: CHEBI:33848 polycyclic arene
"""
from rdkit import Chem

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene based on its SMILES string.
    Polycyclic arenes are hydrocarbons with multiple fused aromatic rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polycyclic arene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find all aromatic (benzene-like) rings
    aromatic_rings = [ring for ring in mol.GetRingInfo().AtomRings() if all(mol.GetAtomWithIdx(atom).GetIsAromatic() for atom in ring)]
    
    # Check if there are at least two aromatic rings
    if len(aromatic_rings) < 2:
        return False, "Less than two aromatic rings found"
    
    # Check for fused rings - a polycyclic structure should have interconnected aromatic rings
    fused_rings = Chem.GetSymmSSSR(mol) >= 2
    if not fused_rings:
        return False, "Aromatic rings are not fused"

    return True, "Contains multiple fused aromatic rings characteristic of polycyclic arenes"

# This function can now be used to classify SMILES strings as polycyclic arenes based on the defined characteristics.
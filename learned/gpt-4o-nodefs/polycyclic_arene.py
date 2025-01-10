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
    
    # Acquire ring information
    ri = mol.GetRingInfo()
    
    # Extract all rings from the molecule, check if aromatic and of decent size
    rings = ri.AtomRings()
    aromatic_rings = [ring for ring in rings if all(mol.GetAtomWithIdx(atom).GetIsAromatic() for atom in ring)]
    
    # Check if there are at least two aromatic rings
    if len(aromatic_rings) < 2:
        return False, "Less than two aromatic rings found"
    
    # Check for fused rings - aromatic rings sharing at least one bond or more
    is_fused = False
    for i, ring1 in enumerate(aromatic_rings):
        for ring2 in aromatic_rings[i+1:]:
            # Check for shared atoms indicating fusion
            if set(ring1).intersection(set(ring2)):
                is_fused = True
                break
        if is_fused:
            break
    
    if not is_fused:
        return False, "Aromatic rings are not fused"

    # Additional check: look for characteristic SMARTS patterns for polycyclic arenes
    # Example SMARTS for a simple polycyclic arene pattern
    polycyclic_smarts = Chem.MolFromSmarts("c1ccccc1c2ccccc2")  # A generic fused bicyclic aromatic ring
    if mol.HasSubstructMatch(polycyclic_smarts):
        return True, "Contains polycyclic arene SMARTS pattern"

    # If none of the aromatic rings are detected with fusion or patterns
    return True, "Contains multiple fused aromatic rings characteristic of polycyclic arenes"
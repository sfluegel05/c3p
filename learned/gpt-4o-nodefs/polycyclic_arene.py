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
    
    # Extract all aromatic rings from the molecule
    aromatic_rings = [ring for ring in ri.AtomRings() if all(mol.GetAtomWithIdx(atom).GetIsAromatic() for atom in ring)]
    
    # Check if there are at least two aromatic rings
    if len(aromatic_rings) < 2:
        return False, "Less than two aromatic rings found"
    
    # Check for fused rings - aromatic rings sharing at least one bond or more
    fused = False
    for i, ring1 in enumerate(aromatic_rings):
        for ring2 in aromatic_rings[i+1:]:
            # Check for shared atoms indicating fusion
            if set(ring1).intersection(set(ring2)):
                fused = True
                break
        if fused:
            break
    
    if not fused:
        return False, "Aromatic rings are not fused"

    return True, "Contains multiple fused aromatic rings characteristic of polycyclic arenes"

# This function can now be used to classify SMILES strings as polycyclic arenes based on the defined characteristics.
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
    
    # Extract all rings from the molecule, verify aromaticity and sizes
    rings = ri.AtomRings()
    aromatic_rings = [ring for ring in rings if Chem.rdMolDescriptors.CalcNumAromaticRings(mol, ring) == len(ring)]
    
    # Check if there are at least two aromatic rings
    if len(aromatic_rings) < 2:
        return False, "Less than two aromatic rings found"
    
    # Check for fused aromatic rings based on shared bonds
    is_fused = False
    fused_rings = []  # store indices of fused rings
    for i, ring1 in enumerate(aromatic_rings):
        for j, ring2 in enumerate(aromatic_rings[i+1:], i+1):
            # Check for shared bonds (edges) not just atoms (nodes)
            if ri.NumBondCrossroads(ring1, ring2) >= 1:
                is_fused = True
                fused_rings.append((ring1, ring2))
        if is_fused:
            break
    
    if not is_fused:
        return False, "Aromatic rings are not fused"
    
    # Additional checks or comprehensive SMARTS patterns could refine the detection here
    # Example: we can define more generalized polycyclic arene patterns
    polycyclic_smarts_patterns = ["c1ccc2cccc3c(=O)ccc4ccccc24c13"]  # examples of polycyclic patterns
    
    for smarts in polycyclic_smarts_patterns:
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            return True, "Matches known polycyclic arene pattern"

    # If none of the aromatic rings are detected with fusion or patterns
    return True, "Contains multiple fused aromatic rings characteristic of polycyclic arenes"
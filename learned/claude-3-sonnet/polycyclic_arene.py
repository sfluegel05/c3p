"""
Classifies: CHEBI:33848 polycyclic arene
"""
"""
Classifies: CHEBI:33510 polycyclic arene
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene based on its SMILES string.
    A polycyclic arene is a polycyclic aromatic hydrocarbon.

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
    
    # Check for aromaticity
    if not mol.GetIsAromatic():
        return False, "The molecule is not aromatic"
    
    # Check for polycyclic structure
    ring_info = mol.GetRingInfo()
    n_rings = len(ring_info.AtomRings())
    if n_rings < 3:
        return False, "The molecule is not polycyclic (fewer than 3 rings)"
    
    # Check for presence of fused rings
    if not ring_info.IsFused():
        return False, "The molecule does not contain fused rings"
    
    # Check for hydrogen deficiency
    hd = rdMolDescriptors.CalcHydrogenDeficiency(mol)
    if hd < 4:
        return False, "Hydrogen deficiency is too low for a polycyclic arene"
    
    # Check for specific substructures
    polycyclic_arene_patterns = [
        Chem.MolFromSmarts("[ar]~[ar]~[ar]~[ar]~[ar]"),
        Chem.MolFromSmarts("[ar]1[ar]2[ar]3[ar]4[ar]5[ar]6[ar]1[ar]2[ar]3[ar]4[ar]5[ar]6"),
        Chem.MolFromSmarts("[ar]1[ar]2[ar]3[ar]4[ar]5[ar]6[ar]7[ar]1[ar]2[ar]3[ar]4[ar]5[ar]6[ar]7")
    ]
    
    has_polycyclic_arene_substructure = any(mol.HasSubstructMatch(pattern) for pattern in polycyclic_arene_patterns)
    
    if not has_polycyclic_arene_substructure:
        return False, "The molecule does not contain the required polycyclic arene substructure"
    
    # If all checks pass, classify as polycyclic arene
    return True, "The molecule is a polycyclic aromatic hydrocarbon"
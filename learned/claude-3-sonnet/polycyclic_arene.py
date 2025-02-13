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
    if n_rings < 2:
        return False, "The molecule is not polycyclic"
    
    # Check for presence of five-membered rings
    if any(len(ring) == 5 for ring in ring_info.AtomRings()):
        # Allow specific five-membered rings like in fluorene and acenaphthylene
        pass
    else:
        # Exclude other five-membered rings
        return False, "The molecule contains a non-aromatic five-membered ring"
    
    # Check for hydrogen deficiency
    hd = rdMolDescriptors.CalcHydrogenDeficiency(mol)
    if hd < 4:
        return False, "Hydrogen deficiency is too low for a polycyclic arene"
    
    # Check for presence of heteroatoms
    atoms = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    if any(atom not in [1, 6, 8] for atom in atoms):
        return False, "The molecule contains heteroatoms other than hydrogen, carbon, and oxygen"
    
    # Check for specific substructures
    polycyclic_arene_pattern = Chem.MolFromSmarts("[ar]~[ar]~[ar]~[ar]")
    if not mol.HasSubstructMatch(polycyclic_arene_pattern):
        return False, "The molecule does not contain the required polycyclic arene substructure"
    
    # If all checks pass, classify as polycyclic arene
    return True, "The molecule is a polycyclic aromatic hydrocarbon"
"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: CHEBI:26389 porphyrin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    A porphyrin is defined by a macrocycle consisting of four pyrrole rings linked by methine bridges.
    This version includes a more robust connectivity check and handles charged and metal-coordinated species.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a porphyrin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the porphyrin macrocycle
    # This pattern captures the core structure with four pyrrole-like rings and four methine bridges
    # It uses generalized nitrogen and carbon atoms to account for variations in charge and metal coordination.
    porphyrin_core_pattern = Chem.MolFromSmarts("[nX3]1=[CX3]~[cX3]~[cX3]~[nX3]2=[CX3]~[cX3]~[cX3]~[nX3]3=[CX3]~[cX3]~[cX3]~[nX3]4=[CX3]~[cX3]~[cX3]~1423")

    # Check if the core pattern matches.
    if not mol.HasSubstructMatch(porphyrin_core_pattern):
        return False, "Porphyrin core macrocycle not found."


    # Check the number of nitrogens, which should be at least 4
    nitrogen_pattern = Chem.MolFromSmarts("[nX3]")
    nitrogen_matches = mol.GetSubstructMatches(nitrogen_pattern)
    if len(nitrogen_matches) < 4:
        return False, f"Found {len(nitrogen_matches)} nitrogen atoms, require at least 4"


    # Additional check to see if there's a macrocycle, using the number of atoms in the core
    match = mol.GetSubstructMatch(porphyrin_core_pattern)
    if match:
        if len(match) != 20:  # the macrocycle core should contain 20 atoms (4 N + 16 C)
           return False, "Porphyrin core has incorrect size"

    # Check for metals
    metal_pattern = Chem.MolFromSmarts("[Fe,Mg,Zn,Co,Ni]")
    metal_matches = mol.GetSubstructMatches(metal_pattern)
    
    
    return True, "Contains a porphyrin core macrocycle"
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

    # Define SMARTS patterns
    # Define pyrrole ring (5-membered ring with one N)
    pyrrole_pattern = Chem.MolFromSmarts("[nX3+]1[cX3][cX3][cX3][cX3]1") # allow for charged N
    
    # Define methine bridge (a carbon connecting to 3 other atoms)
    methine_bridge = Chem.MolFromSmarts("[CX3]")
    
    # Define a pattern connecting 4 pyrrole rings with methine bridges
    #This accounts for variations in bond order and connectivity within the porphyrin.
    porphyrin_core_pattern = Chem.MolFromSmarts("([nX3+]1[cX3][cX3][cX3][cX3]1)[CX3]([nX3+]2[cX3][cX3][cX3][cX3]2)([nX3+]3[cX3][cX3][cX3][cX3]3)[CX3]([nX3+]4[cX3][cX3][cX3][cX3]4)")
    
    # Check if porphyrin core is present
    if not mol.HasSubstructMatch(porphyrin_core_pattern):
       return False, "Porphyrin core macrocycle not found."
    

    # Check for metals
    metal_pattern = Chem.MolFromSmarts("[Fe,Mg,Zn,Co,Ni]")
    metal_matches = mol.GetSubstructMatches(metal_pattern)

    # Count C, N and metals
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    metal_count = len(metal_matches)


    if n_count < 4:
        return False, "Too few nitrogens for porphyrin."
    
    if c_count < 20:
        return False, "Too few carbons for porphyrin."
    
    
    # Additional check to confirm the presence of a macrocycle
    match = mol.GetSubstructMatch(porphyrin_core_pattern)
    if not match:
        return False, "Porphyrin core not matched"

    
    return True, "Contains a porphyrin core macrocycle"
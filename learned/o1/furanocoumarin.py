"""
Classifies: CHEBI:24128 furanocoumarin
"""
"""
Classifies: furanocoumarin
"""

from rdkit import Chem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin consists of a furan ring fused with a coumarin. The fusion
    may occur in different ways, leading to several isomers (linear and angular).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a general SMARTS pattern for the furanocoumarin core
    # Furan ring fused with coumarin (benzopyrone) ring system
    # The tricyclic system can be represented in a general way
    furanocoumarin_smarts = '''
        [$([cR2]1[cR][cR][cR][cR][oR1][cR]1),  # Coumarin ring with oxygen
         $([cR2]1[cR][cR][cR][cR][cR][cR]1)],  # Benzene ring
        $([cR2]1[cR][oR][cR][cR]1)]            # Furan ring
    '''
    furanocoumarin_pattern = Chem.MolFromSmarts(furanocoumarin_smarts)
    if furanocoumarin_pattern is None:
        return False, "Invalid furanocoumarin SMARTS pattern"

    # Find all fused ring systems in the molecule
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()
    fused_rings = []
    for ring in atom_rings:
        if len(ring) >= 5:  # Considering rings of size 5 or more
            fused_rings.append(ring)
    
    # Check for the presence of the tricyclic fused system
    match = mol.HasSubstructMatch(furanocoumarin_pattern)
    if match:
        return True, "Contains furanocoumarin core"
    else:
        return False, "Does not contain furanocoumarin core"
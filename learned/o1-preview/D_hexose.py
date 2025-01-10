"""
Classifies: CHEBI:4194 D-hexose
"""
"""
Classifies: D-hexose
"""

from rdkit import Chem

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose is a hexose (six-carbon monosaccharide) that has D-configuration at position 5.
    This function identifies the C-5 chiral center and checks its configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-hexose, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Assign stereochemistry
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Check for 6 carbons
    num_carbons = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if num_carbons != 6:
        return False, f"Molecule does not have 6 carbons (has {num_carbons})"
    
    # Find all chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
    if len(chiral_centers) != 4:
        return False, f"Molecule does not have 4 chiral centers (has {len(chiral_centers)})"
    
    # Define SMARTS pattern for C-5 chiral center in cyclic form
    # C5 is a chiral carbon connected to:
    # - One oxygen (part of ring)
    # - One carbon (also part of ring)
    # - One carbon (CH2OH group)
    c5_smarts = "[C@H]([CH2][OH])[C]"
    c5_pattern = Chem.MolFromSmarts(c5_smarts)
    matches = mol.GetSubstructMatches(c5_pattern)
    
    if not matches:
        # Try alternative SMARTS pattern for open-chain form
        c5_smarts_linear = "[C@H]([CH2][OH])[CH](O)"
        c5_pattern_linear = Chem.MolFromSmarts(c5_smarts_linear)
        matches = mol.GetSubstructMatches(c5_pattern_linear)
        if not matches:
            return False, "Could not identify C-5 chiral center"
    
    # Check configuration at C-5
    for match in matches:
        c5_idx = match[0]
        atom = mol.GetAtomWithIdx(c5_idx)
        if atom.HasProp('_CIPCode'):
            stereo = atom.GetProp('_CIPCode')
            if stereo == 'R':
                return True, "Molecule is a D-hexose with R configuration at C-5"
            elif stereo == 'S':
                return False, "Molecule is an L-hexose with S configuration at C-5"
            else:
                return False, f"Unknown configuration at C-5: {stereo}"
        else:
            return False, "C-5 chiral center does not have assigned stereochemistry"
    
    return False, "Could not determine configuration at C-5"
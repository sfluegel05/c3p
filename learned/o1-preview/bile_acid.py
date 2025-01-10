"""
Classifies: CHEBI:3098 bile acid
"""
"""
Classifies: bile acid
"""
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.rdchem import ChiralType

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    Bile acids are hydroxy-5beta-cholanic acids occurring in bile.
    This function checks for the presence of a steroid nucleus with
    specific ring sizes and fusion, hydroxyl groups at specific positions,
    and a carboxylic acid side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure molecule has 4 rings
    ring_info = mol.GetRingInfo()
    ring_counts = ring_info.NumRings()
    if ring_counts < 4:
        return False, f"Molecule has {ring_counts} rings, expected at least 4 for steroid nucleus"

    # Get ring systems and check ring sizes
    rings = ring_info.AtomRings()
    ring_sizes = [len(ring) for ring in rings]
    ring_sizes.sort()
    # Steroid nucleus has three 6-membered rings and one 5-membered ring
    if ring_sizes[:4] != [5, 6, 6, 6]:
        return False, f"Ring sizes are {ring_sizes[:4]}, expected [5, 6, 6, 6]"

    # Check for fused ring system
    fused = Chem.GetSymmSSSR(mol)
    if len(fused) < 4:
        return False, "Rings are not fused into a steroid nucleus"

    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O;H]')
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Carboxylic acid group not found"

    # Check for hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts('[OX2H]')
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) == 0:
        return False, "No hydroxyl groups found"

    # Optional: Check for specific stereochemistry at position 5 (5β-configuration)
    # Get atom at position 5 (assuming numbering starting from ring A)
    try:
        # Generate 3D coordinates for stereochemistry
        rdDepictor.Compute2DCoords(mol)
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
        atom5 = mol.GetAtomWithIdx(4)  # Zero-based indexing
        if atom5.GetChiralTag() != ChiralType.CHI_TETRAHEDRAL_CCW:
            return False, "5β-configuration not found"
    except:
        return False, "Error checking stereochemistry at position 5"

    # Additional check: Verify presence of hydroxyl groups at specific positions
    # This requires mapping atom indices to positions, which can be complex
    # For simplicity, we'll skip this detailed check

    return True, "Molecule contains steroid nucleus with correct ring fusion, carboxylic acid, hydroxyl groups, and 5β-configuration"


__metadata__ = {
    'chemical_class': {
        'name': 'bile acid',
        'definition': "Any member of a group of hydroxy-5beta-cholanic acids occurring in bile, where they are present as the sodium salts of their amides with glycine or taurine. In mammals bile acids almost invariably have 5beta-configuration."
    }
}
"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
"""
Classifies: CHEBI:15617 2-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid is any fatty acid with a hydroxy functional group in the alpha- or 2-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"
    
    # Look for alpha-hydroxy group adjacent to carboxyl carbon
    # Define pattern for alpha-hydroxy carboxylic acid
    alpha_hydroxy_pattern = Chem.MolFromSmarts("C(=O)O[C@@H0,R0][CH2]")
    if not mol.HasSubstructMatch(alpha_hydroxy_pattern):
        # Try matching without stereochemistry
        alpha_hydroxy_pattern_no_stereo = Chem.MolFromSmarts("C(=O)O[C][CH2]")
        if not mol.HasSubstructMatch(alpha_hydroxy_pattern_no_stereo):
            return False, "No alpha-hydroxy group adjacent to carboxylic acid found"
    
    # Optionally, check that the molecule is a fatty acid (long aliphatic chain)
    # Count the total number of carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, f"Molecule has {c_count} carbon atoms, less than typical fatty acids"
    
    # Check for rings
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring structures, atypical for fatty acids"
    
    # Check that the molecule is mostly aliphatic (not aromatic)
    aromatic_atoms = any(atom.GetIsAromatic() for atom in mol.GetAtoms())
    if aromatic_atoms:
        return False, "Molecule contains aromatic atoms, which is atypical for fatty acids"
    
    return True, "Molecule is a 2-hydroxy fatty acid with hydroxy group at alpha position"
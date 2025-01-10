"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
"""
Classifies: 3-hydroxy fatty acids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Pattern for 3-hydroxy fatty acids with various possible substitutions at C2
    # Matches both R and S configurations
    beta_hydroxy_pattern = Chem.MolFromSmarts("[OX2H1]-[CX4]-[CX4]-[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(beta_hydroxy_pattern):
        return False, "No hydroxy group at 3-position"

    # Count carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 5:  # Minimum size for 3-hydroxy pentanoic acid
        return False, "Carbon chain too short"
        
    # Count oxygens (should have at least 3: carboxyl group (2) + hydroxy group (1))
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 3:
        return False, "Insufficient oxygen atoms"
        
    # Check for non-allowed atoms (only C, H, O allowed)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6, 8]:
            return False, "Contains non-allowed atoms"
            
    # Check for aromatic systems
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic rings"

    # Get ring information
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        # Check each ring
        ring_atoms = ring_info.AtomRings()
        for ring in ring_atoms:
            if len(ring) != 3:  # Only allow 3-membered rings (cyclopropyl)
                return False, "Contains non-cyclopropyl rings"
    
    # Additional check for the beta-hydroxy pattern positioning
    carboxyl_matches = mol.GetSubstructMatches(carboxylic_pattern)
    beta_hydroxy_matches = mol.GetSubstructMatches(beta_hydroxy_pattern)
    
    if not (carboxyl_matches and beta_hydroxy_matches):
        return False, "Required functional groups not properly positioned"
        
    # Check for proper chain structure
    # Allow branching but ensure main chain exists
    main_chain_pattern = Chem.MolFromSmarts("C-C-C-C")  # Minimum 4 carbons in main chain
    if not mol.HasSubstructMatch(main_chain_pattern):
        return False, "No proper carbon chain found"

    # Verify that hydroxy group is at the 3-position relative to carboxyl group
    # by checking connectivity through SMARTS pattern
    beta_position_pattern = Chem.MolFromSmarts("[OX2H1]-[CX4]-[CX4]-[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(beta_position_pattern):
        return False, "Hydroxy group not properly positioned at 3-position"

    return True, "Molecule contains a 3-hydroxy fatty acid structure"
"""
Classifies: CHEBI:15889 sterol
"""
"""
Classifies: CHEBI:15889 sterol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is a 3-hydroxy steroid whose skeleton is closely related to cholestan-3-ol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic steroid core pattern - four fused rings with flexible bond types
    # Allows for different oxidation states and some variation in structure
    steroid_core = Chem.MolFromSmarts("[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~1")
    
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid ring system found"

    # Check for hydroxyl group (more flexible positioning)
    hydroxyl = Chem.MolFromSmarts("[OX2H1]")
    if not mol.HasSubstructMatches(hydroxyl):
        return False, "No hydroxyl group found"

    # Count carbons (sterols typically have 27-30 carbons but can vary)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbons for a sterol"
    if c_count > 35:
        return False, "Too many carbons for a sterol"

    # Check for angular methyl groups (characteristic of steroids)
    # More flexible pattern that captures different orientations
    angular_methyl = Chem.MolFromSmarts("[CH3][C]([#6])([#6])[#6]")
    methyl_matches = len(mol.GetSubstructMatches(angular_methyl))
    if methyl_matches < 1:
        return False, "Insufficient angular methyl groups"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 700:  # Widened range to catch more variants
        return False, f"Molecular weight {mol_wt:.1f} outside typical sterol range"

    # Check for aliphatic side chain at C-17 (more flexible pattern)
    side_chain = Chem.MolFromSmarts("[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~1[#6]")
    if not mol.HasSubstructMatch(side_chain):
        return False, "No characteristic sterol side chain found"

    # Count rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings"

    # More flexible check for 3-position hydroxyl group
    # Accounts for different possible configurations
    three_pos_oh = Chem.MolFromSmarts("([#6]~1~[#6]~[#6]([OX2H1])~[#6]~[#6]~[#6]~1)|([#6]~1~[#6]~[#6]~[#6]([OX2H1])~[#6]~[#6]~1)")
    if not mol.HasSubstructMatch(three_pos_oh):
        # Try alternative pattern for 3-position
        alt_three_pos_oh = Chem.MolFromSmarts("[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~3~[#6]~[#6]~2~[#6]([OX2H1])~[#6]~1")
        if not mol.HasSubstructMatch(alt_three_pos_oh):
            return False, "No hydroxyl group in characteristic position"

    # Check for reasonable number of oxygens (allowing for more variation)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1 or o_count > 10:
        return False, f"Number of oxygen atoms ({o_count}) outside typical range for sterols"

    return True, "Contains steroid ring system with hydroxyl group and characteristic sterol features"
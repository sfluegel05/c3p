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

    # More flexible steroid core pattern - four fused rings
    # Allows for different bond types and some variation in structure
    steroid_core = Chem.MolFromSmarts("[#6]1[#6][#6][#6]2[#6][#6][#6][#6]3[#6][#6][#6]4[#6][#6][#6][#6][#6]4[#6][#6]3[#6][#6]2[#6]1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid ring system found"

    # Check for hydroxyl group with more flexible positioning
    hydroxyl = Chem.MolFromSmarts("[OH1]")
    if not mol.HasSubstructMatch(hydroxyl):
        return False, "No hydroxyl group found"

    # Count carbons (sterols typically have 27-30 carbons but can vary)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbons for a sterol"
    if c_count > 35:
        return False, "Too many carbons for a sterol"

    # Check for angular methyl groups (characteristic of steroids)
    # More flexible pattern that captures different orientations
    angular_methyls = Chem.MolFromSmarts("[CH3][C]([#6])([#6])[#6]")
    methyl_matches = len(mol.GetSubstructMatches(angular_methyls))
    if methyl_matches < 2:
        return False, f"Found only {methyl_matches} angular methyl groups, need at least 2"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350 or mol_wt > 650:  # Widened range to catch more variants
        return False, f"Molecular weight {mol_wt:.1f} outside typical sterol range"

    # Check for characteristic side chain at C-17
    # More specific pattern that ensures connection to steroid core
    side_chain = Chem.MolFromSmarts("[#6]1[#6][#6][#6]2[#6][#6][#6][#6]3[#6][#6][#6]4[#6][#6][#6][#6][#6]4[#6][#6]3[#6][#6]2[#6]1[#6][#6][#6]")
    if not mol.HasSubstructMatch(side_chain):
        return False, "No characteristic sterol side chain found"

    # Count oxygen atoms (allowing for more variation)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1 or o_count > 8:  # Increased upper limit to allow for more oxidized variants
        return False, f"Number of oxygen atoms ({o_count}) outside typical range for sterols"

    # Check ring count
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings"

    # Additional check for 3-position hydroxyl group
    # More flexible pattern that allows for different configurations
    hydroxy_3_pos = Chem.MolFromSmarts("([#6]1[#6][#6]([OH])[#6][#6][#6]1)|([#6]1[#6][#6][#6]([OH])[#6][#6]1)")
    if not mol.HasSubstructMatch(hydroxy_3_pos):
        return False, "No hydroxyl group in characteristic position"

    return True, "Contains steroid ring system with hydroxyl group and characteristic sterol features"
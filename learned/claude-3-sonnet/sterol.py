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

    # More flexible steroid core pattern that allows for different bond types and modifications
    # This pattern describes the basic 6-6-6-5 tetracyclic system with any type of bonds
    steroid_core = Chem.MolFromSmarts(
        "[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~1"
    )
    
    # Alternative core pattern that's more flexible with ring fusion points
    alt_core = Chem.MolFromSmarts(
        "[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~1"
    )

    if not (mol.HasSubstructMatch(steroid_core) or mol.HasSubstructMatch(alt_core)):
        return False, "No steroid ring system found"

    # Check for any hydroxyl group
    hydroxyl = Chem.MolFromSmarts("[OX2H1]")
    if not mol.HasSubstructMatches(hydroxyl):
        return False, "No hydroxyl group found"

    # Count carbons (sterols typically have 27-30 carbons but can vary)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbons for a sterol"
    if c_count > 40:  # Increased upper limit to account for more complex sterols
        return False, "Too many carbons for a sterol"

    # Check for angular methyl groups (characteristic of steroids)
    # More flexible pattern that captures different orientations
    angular_methyl = Chem.MolFromSmarts("[CH3][C]([#6])([#6])[#6]")
    methyl_matches = len(mol.GetSubstructMatches(angular_methyl))
    if methyl_matches < 1:
        return False, "Insufficient angular methyl groups"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 800:  # Widened range further to catch more variants
        return False, f"Molecular weight {mol_wt:.1f} outside typical sterol range"

    # Count rings (sterols should have at least 4 rings)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings"

    # More flexible check for hydroxyl group position
    # This pattern looks for a hydroxyl connected to any carbon in the ring system
    ring_hydroxyl = Chem.MolFromSmarts("([#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1)[#6][OH]")
    alt_ring_hydroxyl = Chem.MolFromSmarts("([#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1)~[#6]~[OH]")
    
    if not (mol.HasSubstructMatch(ring_hydroxyl) or mol.HasSubstructMatch(alt_ring_hydroxyl)):
        return False, "No hydroxyl group in characteristic position"

    # Count oxygens (allowing for more variation)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1 or o_count > 12:  # Increased upper limit to account for more oxidized sterols
        return False, f"Number of oxygen atoms ({o_count}) outside typical range for sterols"

    # Check for reasonable degree of unsaturation
    double_bonds = Chem.MolFromSmarts("C=C")
    n_double_bonds = len(mol.GetSubstructMatches(double_bonds))
    if n_double_bonds > 8:  # Sterols typically don't have too many double bonds
        return False, "Too many double bonds for a sterol"

    return True, "Contains steroid ring system with hydroxyl group and characteristic sterol features"
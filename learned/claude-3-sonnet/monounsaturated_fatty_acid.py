"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
"""
Classifies: monounsaturated fatty acids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid based on its SMILES string.
    Must have exactly one double/triple bond in the main chain and a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    if carboxylic_pattern is None or carboxylate_pattern is None:
        return None, None
    
    has_acid = mol.HasSubstructMatch(carboxylic_pattern)
    has_carboxylate = mol.HasSubstructMatch(carboxylate_pattern)
    
    if not (has_acid or has_carboxylate):
        return False, "No carboxylic acid or carboxylate group found"

    # Count carbons in molecule
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Carbon chain too short for fatty acid"

    # Count different types of unsaturations
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    cyclopropene_pattern = Chem.MolFromSmarts("C1=CC1")
    
    if any(pattern is None for pattern in [double_bond_pattern, triple_bond_pattern, cyclopropene_pattern]):
        return None, None
    
    double_bonds = len(mol.GetSubstructMatches(double_bond_pattern))
    triple_bonds = len(mol.GetSubstructMatches(triple_bond_pattern))
    cyclopropene_rings = len(mol.GetSubstructMatches(cyclopropene_pattern))
    
    total_unsaturations = double_bonds + triple_bonds + cyclopropene_rings
    
    if total_unsaturations == 0:
        return False, "No unsaturations found"
    elif total_unsaturations > 1:
        # Special case for cyclopropene-containing fatty acids
        if cyclopropene_rings == 1 and double_bonds == 0 and triple_bonds == 0:
            return True, "Contains one unsaturation in cyclopropene ring"
        return False, f"Too many unsaturations: {total_unsaturations}"

    # Check for ring structures
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    
    # Allow 3-membered rings (cyclopropene, cyclopropane, epoxide)
    for size in ring_sizes:
        if size > 3 and size != 5:  # Allow 3-membered rings and 5-membered rings
            return False, f"Contains unexpected ring of size {size}"

    # Check for reasonable chain length
    chain_pattern = Chem.MolFromSmarts("C(=O)[OH,O-]~C~C~C")
    if chain_pattern is None:
        return None, None
    
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No proper fatty acid chain found"

    # Additional check for molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 70:  # Minimum weight for but-3-enoic acid
        return False, "Molecule too small for fatty acid"
    if mol_wt > 800:  # Allow for some derivatives and larger fatty acids
        return False, "Molecule too large for fatty acid"

    # Check for common substituents that are allowed
    allowed_substituents = [
        ("O", "hydroxy"),
        ("OC(=O)", "ester"),
        ("[O-]", "alkoxide")
    ]
    
    for pattern, _ in allowed_substituents:
        subst_pattern = Chem.MolFromSmarts(pattern)
        if subst_pattern is None:
            return None, None

    # Check for forbidden groups that would make it not a fatty acid
    forbidden_groups = [
        ("N", "nitrogen"),
        ("S", "sulfur"),
        ("B", "boron"),
        ("P", "phosphorus"),
        ("Si", "silicon"),
        ("F", "fluorine"),
        ("Cl", "chlorine"),
        ("Br", "bromine"),
        ("I", "iodine")
    ]
    
    for pattern, group_name in forbidden_groups:
        group_pattern = Chem.MolFromSmarts(pattern)
        if group_pattern is None:
            return None, None
        if mol.HasSubstructMatch(group_pattern):
            return False, f"Contains forbidden group: {group_name}"

    return True, "Contains one unsaturation and appropriate carboxylic acid group in fatty acid structure"
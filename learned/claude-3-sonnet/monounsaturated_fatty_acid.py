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

    # Check for carboxylic acid or carboxylate group
    carboxylic_pattern = Chem.MolFromSmarts("[$([CX3](=[OX1])[OX2H1]),$([CX3](=[OX1])[OX2-])]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid or carboxylate group found"

    # Count carbons in molecule
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Carbon chain too short for fatty acid"

    # Count double and triple bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    cyclopropene_pattern = Chem.MolFromSmarts("C1=CC1") # For sterculic acid type structures
    
    double_bonds = len(mol.GetSubstructMatches(double_bond_pattern))
    triple_bonds = len(mol.GetSubstructMatches(triple_bond_pattern))
    cyclopropene_rings = len(mol.GetSubstructMatches(cyclopropene_pattern))
    
    # Count total unsaturations (cyclopropene counts as one unsaturation)
    total_unsaturations = double_bonds + triple_bonds + cyclopropene_rings
    
    if total_unsaturations == 0:
        return False, "No unsaturations found"
    elif total_unsaturations > 1:
        # Special case: check if one unsaturation is in a cyclopropene ring
        if cyclopropene_rings == 1 and double_bonds == 0 and triple_bonds == 0:
            return True, "Contains one unsaturation in cyclopropene ring (sterculic acid type)"
        return False, f"Too many unsaturations: {total_unsaturations}"

    # Check for conjugated systems
    conjugated_pattern = Chem.MolFromSmarts("C=C-C=C")
    if mol.HasSubstructMatch(conjugated_pattern):
        return False, "Contains conjugated double bonds"

    # Check for forbidden groups
    forbidden_groups = [
        ("C(=O)N", "amide"),
        ("S(=O)(=O)", "sulfone/sulfoxide"),
        ("B", "boron"),
        ("P", "phosphorus"),
        ("Si", "silicon"),
    ]
    
    for pattern, group_name in forbidden_groups:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return False, f"Contains forbidden group: {group_name}"

    # Check for reasonable chain length
    # Look for carbon chain connected to carboxylic group
    chain_pattern = Chem.MolFromSmarts("C(=O)O[H,-]~C~C~C")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No proper fatty acid chain found"

    # Allow common substituents and check ring structures
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    
    # Allow small rings (cyclopropane, epoxide, cyclopropene) but reject larger ones
    # except for specific cases like sterculic acid
    for size in ring_sizes:
        if size > 3 and size != 5:  # Allow 3-membered rings and 5-membered rings
            return False, f"Contains unexpected ring of size {size}"

    # Additional check for molecular weight to ensure reasonable size
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 70:  # Minimum weight for but-3-enoic acid
        return False, "Molecule too small for fatty acid"
    if mol_wt > 400:  # Maximum weight for typical monounsaturated fatty acids
        # Check if it's a derivative with allowed substituents
        if not any(mol.HasSubstructMatch(Chem.MolFromSmarts(p)) for p in ["[OH]", "[O-]", "OC(=O)"]):
            return False, "Molecule too large for typical fatty acid"

    return True, "Contains one unsaturation and appropriate carboxylic acid group in fatty acid structure"
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
    Must have exactly one double/triple bond in the main chain and a carboxylic acid group
    or its derivatives.

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

    # Count carbons in main chain
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Carbon chain too short for fatty acid"

    # Count unsaturations in carbon chain (excluding carboxylic group)
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    
    double_bonds = len(mol.GetSubstructMatches(double_bond_pattern))
    triple_bonds = len(mol.GetSubstructMatches(triple_bond_pattern))
    
    total_unsaturations = double_bonds + triple_bonds
    
    if total_unsaturations == 0:
        return False, "No unsaturations found in carbon chain"
    elif total_unsaturations > 1:
        # Check if the additional unsaturations are part of other functional groups
        conjugated_pattern = Chem.MolFromSmarts("C=C-C=C")
        if mol.HasSubstructMatch(conjugated_pattern):
            return False, "Contains conjugated double bonds"
            
        # Check for unsaturations in substituents vs main chain
        chain_pattern = Chem.MolFromSmarts("C(=O)O[H,\-]~C~C~C~C")
        if not mol.HasSubstructMatch(chain_pattern):
            return False, "Multiple unsaturations in main chain"
    
    # Allow common substituents
    allowed_substituents = [
        Chem.MolFromSmarts("[OX2H1]"), # hydroxyl
        Chem.MolFromSmarts("[OX2][CX4]"), # alkoxy
        Chem.MolFromSmarts("[NX3H2]"), # amino
        Chem.MolFromSmarts("[F,Cl,Br,I]") # halogens
    ]
    
    # Count non-allowed functional groups
    forbidden_groups = [
        Chem.MolFromSmarts("C(=O)N"), # amide
        Chem.MolFromSmarts("S(=O)(=O)"), # sulfone/sulfoxide
        Chem.MolFromSmarts("B"), # boron
        Chem.MolFromSmarts("P"), # phosphorus
        Chem.MolFromSmarts("Si"), # silicon
    ]
    
    for pattern in forbidden_groups:
        if pattern and mol.HasSubstructMatch(pattern):
            return False, "Contains non-fatty acid functional groups"

    # Check for reasonable chain length and structure
    chain_length = max(len(path) for path in rdMolDescriptors.FindAllPathsOfLengthN(mol, 8))
    if chain_length < 5:
        return False, "Carbon chain too short or disconnected"
        
    # Allow small rings (cyclopropane, epoxide) but reject larger ones
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    if any(size > 4 for size in ring_sizes):
        # Special check for common fatty acid rings
        if any(size == 5 for size in ring_sizes) and c_count <= 20:
            return True, "Contains cyclopentyl ring structure common in some fatty acids"
        return False, "Contains large ring structure"

    return True, "Contains appropriate unsaturation and carboxylic acid group in fatty acid structure"
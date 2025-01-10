"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
"""
Classifies: trienoic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid is any polyunsaturated fatty acid that contains three double bonds.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a trienoic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that molecule contains only C, H, and O atoms
    allowed_atoms = {1, 6, 8}  # H, C, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, f"Molecule contains disallowed atom {atom.GetSymbol()}"
    
    # Check for exactly one carboxylic acid group (C(=O)O[H])
    carboxylic_acid = Chem.MolFromSmarts("C(=O)[O;H1]")
    matches = mol.GetSubstructMatches(carboxylic_acid)
    if len(matches) != 1:
        return False, f"Found {len(matches)} carboxylic acid groups, expected exactly 1"
    
    # Count total number of carbons
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 12:
        return False, f"Molecule contains {num_carbons} carbons, expected at least 12"
    
    # Count number of aliphatic double bonds
    num_double_bonds = rdMolDescriptors.CalcNumAliphaticDoubleBonds(mol)
    # Count number of aliphatic rings
    num_rings = rdMolDescriptors.CalcNumAliphaticRings(mol)
    
    total_unsaturation = num_double_bonds + num_rings
    
    if total_unsaturation < 3:
        return False, f"Contains {num_double_bonds} aliphatic double bonds and {num_rings} aliphatic rings, total unsaturation {total_unsaturation}, expected at least 3"
    
    return True, "Molecule is a trienoic fatty acid: contains exactly one carboxylic acid group, at least 12 carbons, only C, H, O atoms, and total unsaturation (double bonds + rings) of at least 3"

__metadata__ = {
    'chemical_class': {
        'name': 'trienoic fatty acid',
        'definition': 'Any polyunsaturated fatty acid that contains three double bonds.',
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    }
}
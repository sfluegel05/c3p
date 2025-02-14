"""
Classifies: CHEBI:17761 ceramide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a ceramide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Key patterns to identify
    sphingoid_base_pattern = Chem.MolFromSmarts("[NX3][C@@H](CO)[C@H](O)[C@@H](O)C")  # General sphingosine base
    amide_fatty_acid_pattern = Chem.MolFromSmarts("C(=O)[NX3][C@H](C)C")  # Amide linkage with aliphatic extension
    hydroxyl_group_pattern = Chem.MolFromSmarts("[O]")
    
    # Match patterns
    if not mol.HasSubstructMatch(sphingoid_base_pattern):
        return False, "No sphingoid base pattern found"
    
    if not mol.HasSubstructMatch(amide_fatty_acid_pattern):
        return False, "No amide-linked fatty acid found"
    
    # Ensure presence of at least one hydroxyl group
    hydroxyl_groups = mol.GetSubstructMatches(hydroxyl_group_pattern)
    if len(hydroxyl_groups) < 1:
        return False, "No hydroxyl group found"
    
    # Get carbon chain lengths
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_heteroatoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() != 6 and atom.GetAtomicNum() != 1)
    
    if c_count < 18 or c_count > 30:
        return False, f"Total carbon count {c_count} out of typical ceramide range 18-30"

    # Determine molecular weight to approximate fatty acid chain length
    mol_wt = Descriptors.MolWt(mol)
    if mol_wt < 400 or mol_wt > 1000:
        return False, "Molecular weight out of typical ceramide range"

    return True, "Contains sphingoid base with amide-linked fatty acid chain meeting ceramide criteria"
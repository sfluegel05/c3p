"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
"""
Classifies: CHEBI:17895 dicarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid has exactly two carboxy (-COOH) groups and no anhydride bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dicarboxylic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carboxylic acid groups using RDKit's built-in method
    num_carboxy = rdMolDescriptors.CalcNumCarboxylicAcids(mol)
    
    if num_carboxy != 2:
        return False, f"Found {num_carboxy} carboxylic acid groups, need exactly 2"

    # Check for anhydride pattern (two carbonyls connected via oxygen)
    anhydride_pattern = Chem.MolFromSmarts("[CX3](=O)-O-[CX3](=O)")
    if mol.HasSubstructMatch(anhydride_pattern):
        return False, "Contains anhydride linkage between carboxy groups"

    # Verify groups are distinct (optional check for conjugated systems)
    carboxy_atoms = set()
    for group in mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3](=O)[OH]")):
        carboxy_atoms.update(group)
    
    # Check we have exactly 4 oxygens in carboxy groups (2 groups * 2 oxygens each)
    o_in_groups = sum(1 for atom in mol.GetAtoms() 
                      if atom.GetAtomicNum() == 8 and atom.GetIdx() in carboxy_atoms)
    
    if o_in_groups != 4:
        return False, "Oxygen count mismatch in carboxy groups"

    return True, "Contains two distinct carboxylic acid groups"
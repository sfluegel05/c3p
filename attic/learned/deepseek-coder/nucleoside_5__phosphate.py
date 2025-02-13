"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
"""
Classifies: CHEBI:16761 nucleoside 5'-phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    A nucleoside 5'-phosphate consists of a ribose or deoxyribose sugar with a purine or pyrimidine base
    and a phosphate group attached to the 5' carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside 5'-phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define patterns for key components
    # Ribose or deoxyribose sugar ring pattern
    sugar_pattern = Chem.MolFromSmarts("[C@H]1O[C@H](COP(=O)(O)O)[C@@H](O)[C@H]1O")
    
    # Purine base pattern (adenine, guanine)
    purine_pattern = Chem.MolFromSmarts("n1cnc2ncnc12")
    
    # Pyrimidine base pattern (cytosine, thymine, uracil)
    pyrimidine_pattern = Chem.MolFromSmarts("n1ccc(=O)[nH]c1=O")
    
    # Check for sugar ring with phosphate group
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No ribose/deoxyribose sugar with 5'-phosphate found"
    
    # Check for purine or pyrimidine base
    has_purine = mol.HasSubstructMatch(purine_pattern)
    has_pyrimidine = mol.HasSubstructMatch(pyrimidine_pattern)
    
    if not (has_purine or has_pyrimidine):
        return False, "No purine or pyrimidine base found"
    
    # Count phosphate groups
    phosphate_pattern = Chem.MolFromSmarts("[PX4](=O)(O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 1:
        return False, "No phosphate group found"
    
    # Check molecular weight - nucleoside 5'-phosphates typically >300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for nucleoside 5'-phosphate"
    
    # Count carbons, oxygens, and nitrogens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if c_count < 5:
        return False, "Too few carbons for nucleoside 5'-phosphate"
    if o_count < 5:
        return False, "Too few oxygens for nucleoside 5'-phosphate"
    if n_count < 1:
        return False, "No nitrogen found in base"
    
    return True, "Contains ribose/deoxyribose sugar with purine/pyrimidine base and 5'-phosphate group"
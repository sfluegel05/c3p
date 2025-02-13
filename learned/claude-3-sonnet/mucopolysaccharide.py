"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
"""
Classifies: CHEBI:26611 mucopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide based on its SMILES string.
    A mucopolysaccharide is a polysaccharide composed of alternating uronic acids and glycosamines,
    and commonly partially esterified with sulfuric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mucopolysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for uronic acid and glycosamine units
    uronic_acid_pattern = Chem.MolFromSmarts("[C&X3](=[O&X1])([O&X2])[O&X2]")
    glycosamine_pattern = Chem.MolFromSmarts("[N&X3]([C&X4])[C&X4]([O&X2])[C&X4]([O&X2])[C&X4]([O&X2])[C&X4]")
    
    # Check for alternating uronic acid-glycosamine pattern
    uronic_acid_matches = mol.GetSubstructMatches(uronic_acid_pattern)
    glycosamine_matches = mol.GetSubstructMatches(glycosamine_pattern)
    
    if not uronic_acid_matches or not glycosamine_matches:
        return False, "Missing uronic acid or glycosamine units"
    
    alternating_pattern = False
    for uronic_acid_idx in uronic_acid_matches:
        for glycosamine_idx in glycosamine_matches:
            if abs(uronic_acid_idx - glycosamine_idx) == 1:
                alternating_pattern = True
                break
        if alternating_pattern:
            break
    
    if not alternating_pattern:
        return False, "Uronic acid and glycosamine units not in alternating pattern"
    
    # Look for sulfate groups or sulfur atoms
    has_sulfate = any(atom.GetAtomicNum() == 16 for atom in mol.GetAtoms())
    
    # Check molecular weight and atom counts
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for mucopolysaccharide"
    
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 10 or o_count < 5:
        return False, "Insufficient carbon and oxygen atoms for polysaccharide"
    
    # Classify as mucopolysaccharide if alternating pattern and sulfation present
    if alternating_pattern and has_sulfate:
        return True, "Contains alternating uronic acid and glycosamine units, partially esterified with sulfate groups"
    else:
        return False, "Does not match mucopolysaccharide structural features"
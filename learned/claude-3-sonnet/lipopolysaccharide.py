"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
"""
Classifies: CHEBI:16852 lipopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    Lipopolysaccharides are complex molecules containing both sugar and lipid components.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopolysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate basic properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:  # Lipopolysaccharides are typically large molecules
        return False, "Molecular weight too low for lipopolysaccharide"

    # Count key atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if o_count < 6:  # Need multiple oxygen atoms for glycosidic bonds and hydroxyl groups
        return False, "Too few oxygen atoms for lipopolysaccharide"
    
    if c_count < 12:  # Need carbon chains for both sugar and lipid components
        return False, "Too few carbon atoms for lipopolysaccharide"

    # Look for sugar ring patterns (pyranose or furanose)
    sugar_pattern = Chem.MolFromSmarts("[CR1][OR1][CR1][CR1][CR1][CR1]")  # 6-membered sugar ring
    sugar_pattern2 = Chem.MolFromSmarts("[CR1][OR1][CR1][CR1][CR1]")  # 5-membered sugar ring
    sugar_matches = len(mol.GetSubstructMatches(sugar_pattern)) + len(mol.GetSubstructMatches(sugar_pattern2))
    
    if sugar_matches < 1:
        return False, "No sugar rings found"

    # Look for hydroxyl groups (characteristic of sugars)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    if hydroxyl_matches < 3:
        return False, "Too few hydroxyl groups for lipopolysaccharide"

    # Look for lipid characteristics (long carbon chains)
    alkyl_chain = Chem.MolFromSmarts("[CH2][CH2][CH2][CH2]")
    alkyl_matches = len(mol.GetSubstructMatches(alkyl_chain))
    
    if alkyl_matches < 1:
        return False, "No long alkyl chains found"

    # Look for glycosidic bonds or ester linkages
    glycosidic_pattern = Chem.MolFromSmarts("[CX4][OX2][CX4]")
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2]")
    linkage_matches = len(mol.GetSubstructMatches(glycosidic_pattern)) + len(mol.GetSubstructMatches(ester_pattern))
    
    if linkage_matches < 2:
        return False, "Insufficient glycosidic/ester linkages"

    # Count rotatable bonds to verify structural complexity
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Structure too rigid for lipopolysaccharide"

    # Check for minimum structural complexity
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 1:
        return False, "No rings found"

    return True, "Contains sugar units, lipid components, and appropriate linkages characteristic of lipopolysaccharides"
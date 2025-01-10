"""
Classifies: CHEBI:25106 macrolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.
    A macrolide is a macrocyclic lactone with a ring of twelve or more members derived from a polyketide.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a macrolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Look for ester group in ring (-O-C(=O)-)
    lactone_pattern = Chem.MolFromSmarts("[O;R][C;R](=[O;!R])")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone group found"
    
    # Find the macrocyclic ring containing the lactone
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings():
        return False, "No rings found"
    
    # Get ring sizes containing the lactone
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    max_ring_size = 0
    
    for match in lactone_matches:
        # Get all rings containing these atoms
        for ring in ring_info.AtomRings():
            if match[0] in ring and match[1] in ring:  # if lactone atoms are in ring
                max_ring_size = max(max_ring_size, len(ring))
    
    if max_ring_size < 12:
        return False, f"Largest ring containing lactone has {max_ring_size} members, needs 12+"
        
    # Check for typical polyketide-derived features
    
    # Count carbonyls (ketones and esters)
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])")
    carbonyl_count = len(mol.GetSubstructMatches(carbonyl_pattern))
    
    # Count hydroxyls
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    # Count methyl groups
    methyl_pattern = Chem.MolFromSmarts("[CH3]")
    methyl_count = len(mol.GetSubstructMatches(methyl_pattern))
    
    # Typical macrolides should have multiple of these features
    feature_count = (carbonyl_count > 1) + (hydroxyl_count > 0) + (methyl_count > 0)
    if feature_count < 2:
        return False, "Insufficient polyketide-like features"
    
    # Check molecular weight - most macrolides are fairly large
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:  # arbitrary cutoff but most macrolides are larger
        return False, "Molecular weight too low for typical macrolide"
        
    # Count sp3 carbons - macrolides typically have many
    sp3_c = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum()==6 and atom.GetHybridization()==Chem.HybridizationType.SP3)
    if sp3_c < 5:
        return False, "Too few sp3 carbons for typical macrolide"
    
    return True, f"Contains {max_ring_size}-membered macrocyclic lactone with polyketide features"
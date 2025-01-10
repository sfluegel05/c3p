"""
Classifies: CHEBI:51963 hopanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_hopanoid(smiles: str):
    """
    Determines if a molecule is a hopanoid based on its SMILES string.
    Hopanoids are pentacyclic triterpenoids based on the hopane skeleton.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a hopanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Hopane core SMARTS patterns - using multiple patterns to catch different variations
    # Pattern 1: Basic hopane skeleton with flexible ring connections
    hopane_core1 = Chem.MolFromSmarts("[C,c]1[C,c]2[C,c][C,c][C,c]3[C,c]([C,c]2[C,c][C,c][C,c]1)[C,c][C,c][C,c]4[C,c]3[C,c][C,c][C,c]4")
    
    # Pattern 2: Alternative hopane skeleton pattern focusing on ring fusion points
    hopane_core2 = Chem.MolFromSmarts("[C,c]1[C,c][C,c][C,c]2[C,c]3[C,c][C,c][C,c]4[C,c]([C,c]3[C,c][C,c]2[C,c]1)[C,c][C,c][C,c][C,c]4")

    if not (mol.HasSubstructMatch(hopane_core1) or mol.HasSubstructMatch(hopane_core2)):
        return False, "No hopane pentacyclic core found"

    # Count rings - hopanoids should have at least 5 rings
    ri = mol.GetRingInfo()
    if ri.NumRings() < 5:
        return False, "Insufficient number of rings"

    # Check for characteristic angular methyl groups at key positions
    angular_methyl = Chem.MolFromSmarts("[C]1[C]([CH3])[C][C]2")
    if len(mol.GetSubstructMatches(angular_methyl)) < 2:
        return False, "Missing characteristic angular methyl groups"

    # Verify basic molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350:  # Remove upper limit to allow for modified hopanoids
        return False, f"Molecular weight ({mol_wt:.1f}) too low for hopanoid"

    # Additional structural features to differentiate from other triterpenes
    
    # Check for characteristic E ring pattern (6-membered)
    e_ring = Chem.MolFromSmarts("[C]1[C][C][C][C]([C]1)([C,CH3])[C,CH3]")
    if not mol.HasSubstructMatch(e_ring):
        return False, "Missing characteristic E ring pattern"

    # Check for gem-dimethyl group often present in hopanoids
    gem_dimethyl = Chem.MolFromSmarts("[C]([CH3])([CH3])[C,CH2,CH]")
    if not mol.HasSubstructMatch(gem_dimethyl):
        return False, "Missing characteristic gem-dimethyl group"

    # Look for characteristic hopanoid side chain patterns
    side_chain_patterns = [
        Chem.MolFromSmarts("C(C)(C)O"),  # Common hydroxylated side chain
        Chem.MolFromSmarts("C(C)(C)CC"),  # Extended side chain
        Chem.MolFromSmarts("C(=C)(C)C"),  # Unsaturated side chain
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in side_chain_patterns):
        return False, "Missing characteristic side chain patterns"

    return True, "Contains hopane pentacyclic core with characteristic features"
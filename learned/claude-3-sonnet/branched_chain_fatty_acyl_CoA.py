"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
"""
Classifies: branched-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a branched-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA moiety components
    
    # Look for adenine pattern
    adenine_pattern = Chem.MolFromSmarts("c1nc2c(n1)c(N)ncn2")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "No CoA moiety found (missing adenine)"
    
    # Look for thioester linkage (-C(=O)S-)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"
        
    # Look for phosphate groups characteristic of CoA
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 2:
        return False, "Missing characteristic CoA phosphate groups"

    # Check for various types of branching patterns
    
    # Standard methyl branch
    methyl_branch = Chem.MolFromSmarts("[CH3][CH1](-[#6])-[#6]")
    
    # Isopropyl group (as in leucine derivatives)
    isopropyl = Chem.MolFromSmarts("[CH3][CH]([CH3])[CH2]")
    
    # Secondary methyl branch (as in iso-compounds)
    iso_branch = Chem.MolFromSmarts("[CH3][CH](-[#6])-[#6]")
    
    # Tertiary carbon branch
    tert_branch = Chem.MolFromSmarts("[#6][C](-[#6])(-[#6])-[#6]")
    
    # Multiple methyl groups on same carbon
    geminal_dimethyl = Chem.MolFromSmarts("[CH3][C]([CH3])(-[#6])-[#6]")
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in 
              [methyl_branch, isopropyl, iso_branch, tert_branch, geminal_dimethyl]):
        return False, "No branching patterns found in carbon chain"

    # Verify chain length - need at least 4 carbons in main chain
    # (excluding CoA part) connected to the thioester
    chain_pattern = Chem.MolFromSmarts("[#6]~[#6]~[#6]~[#6]~[CX3](=O)[SX2]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Carbon chain too short for branched-chain fatty acid"

    # Count carbons to ensure reasonable size
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 23:  # CoA itself contributes about 23 carbons
        return False, "Total carbon count too low for branched-chain fatty acyl-CoA"

    # Optional: Check for common modifications
    hydroxy_pattern = Chem.MolFromSmarts("[#6][CH1]([OH1])[#6]")
    oxo_pattern = Chem.MolFromSmarts("[#6]C(=O)[#6]")
    
    modifications = []
    if mol.HasSubstructMatch(hydroxy_pattern):
        modifications.append("hydroxy-modified")
    if mol.HasSubstructMatch(oxo_pattern):
        modifications.append("oxo-modified")
        
    base_reason = "Contains CoA moiety, thioester linkage, and branched fatty acid chain"
    if modifications:
        return True, f"{base_reason} ({', '.join(modifications)})"
    return True, base_reason
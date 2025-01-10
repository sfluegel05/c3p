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
    adenine_pattern = Chem.MolFromSmarts("c1nc2c(n1)c(N)ncn2")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "No CoA moiety found (missing adenine)"
    
    # Look for thioester linkage (-C(=O)S-)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    if not mol.HasSubstructMatches(thioester_pattern):
        return False, "No thioester linkage found"
        
    # Look for phosphate groups characteristic of CoA
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 2:
        return False, "Missing characteristic CoA phosphate groups"

    # Check for cyclic structures in the fatty acid part
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 1:  # Allow 1 ring for adenine
        return False, "Contains additional ring structures - not a simple branched chain"

    # Define specific branching patterns for branched-chain fatty acids
    branching_patterns = [
        # Isopropyl terminus (as in leucine derivatives)
        Chem.MolFromSmarts("[CH3][CH]([CH3])[CH2][CH2]C(=O)[S]"),
        
        # Isobutyryl (from valine)
        Chem.MolFromSmarts("[CH3][CH]([CH3])C(=O)[S]"),
        
        # 2-methylbutyryl (from isoleucine)
        Chem.MolFromSmarts("[CH3][CH2][CH]([CH3])C(=O)[S]"),
        
        # Other methyl branches along chain
        Chem.MolFromSmarts("[CH3][CH]([CH2][CH2])C(=O)[S]"),
        
        # Longer iso-branches
        Chem.MolFromSmarts("[CH3][CH]([CH2][CH2][CH2])C(=O)[S]"),
        
        # Multiple methyl branches
        Chem.MolFromSmarts("[CH3][CH]([CH2][CH2])[CH]([CH3])")
    ]

    found_branching = False
    for pattern in branching_patterns:
        if mol.HasSubstructMatch(pattern):
            found_branching = True
            break
            
    if not found_branching:
        return False, "No characteristic branched-chain fatty acid pattern found"

    # Optional: Check for common modifications that don't affect classification
    modifications = []
    
    # Check for hydroxy groups
    hydroxy_pattern = Chem.MolFromSmarts("[CH1]([OH1])")
    if mol.HasSubstructMatch(hydroxy_pattern):
        modifications.append("hydroxy-modified")
        
    # Check for additional oxo groups (besides thioester)
    oxo_pattern = Chem.MolFromSmarts("[CH2]C(=O)[CH2]")
    if mol.HasSubstructMatch(oxo_pattern):
        modifications.append("oxo-modified")

    base_reason = "Contains CoA moiety, thioester linkage, and characteristic branched-chain fatty acid pattern"
    if modifications:
        return True, f"{base_reason} ({', '.join(modifications)})"
    return True, base_reason
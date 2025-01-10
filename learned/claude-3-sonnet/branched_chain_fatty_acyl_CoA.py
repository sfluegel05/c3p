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

    # Check for CoA moiety
    # Look for adenine pattern
    adenine_pattern = Chem.MolFromSmarts("c1nc2c(n1)c(N)ncn2")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "No CoA moiety found (missing adenine)"
    
    # Look for thioester linkage (-C(=O)S-)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    if not mol.HasSubstructMatches(thioester_pattern):
        return False, "No thioester linkage found"
        
    # Look for phosphate groups characteristic of CoA
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]")
    if len(mol.GetSubstructMatches(phosphate_pattern)) < 2:
        return False, "Missing characteristic CoA phosphate groups"

    # Check for branched chain:
    # Look for carbons with more than 2 carbon neighbors (branch points)
    branch_pattern = Chem.MolFromSmarts("[CH1,CH0](-[#6])(-[#6])(-[#6])")
    if not mol.HasSubstructMatch(branch_pattern):
        return False, "No branching found in carbon chain"

    # Look for methyl branches specifically
    methyl_branch_pattern = Chem.MolFromSmarts("[CH3][CH1]")
    if not mol.HasSubstructMatch(methyl_branch_pattern):
        iso_pattern = Chem.MolFromSmarts("[CH3][CH]([CH3])[CH2]")
        if not mol.HasSubstructMatch(iso_pattern):
            return False, "No methyl branches or isopropyl groups found"

    # Check for minimum chain length (typically these are fatty acids)
    carbon_chain = Chem.MolFromSmarts("[#6]~[#6]~[#6]~[#6]")
    if not mol.HasSubstructMatch(carbon_chain):
        return False, "Carbon chain too short"

    # Count carbons to verify it's a fatty acid derivative
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:  # CoA itself contributes many carbons
        return False, "Total carbon count too low for branched-chain fatty acyl-CoA"

    return True, "Contains CoA moiety, thioester linkage, and branched fatty acid chain"
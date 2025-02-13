"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
"""
Classifies: medium-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA(4-) based on its SMILES string.
    Medium-chain fatty acids typically have 6-12 carbons in their chain.
    The (4-) species has exactly 4 negative charges.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for CoA phosphate groups - using more flexible pattern
    phosphate_pattern = Chem.MolFromSmarts("[PX4](=[OX1])([OX2H,OX1-])[OX2H,OX1-]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 3:
        return False, f"Found {len(phosphate_matches)} phosphate groups, need at least 3"
        
    # Count negative charges
    explicit_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if explicit_charge != -4:
        return False, f"Total charge is {explicit_charge}, need -4"
        
    # Check for adenine
    adenine_pattern = Chem.MolFromSmarts("c1nc(c2c(n1)ncn2)N")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Missing adenine moiety"
        
    # Check for thioester bond
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Missing thioester bond"
        
    # Check for pantetheine arm
    pantetheine_pattern = Chem.MolFromSmarts("NCCC(=O)NCCS")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Missing pantetheine arm"

    # Count carbons in fatty acid chain
    # CoA has 23 carbons, so subtract to get fatty acid chain length
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    fatty_acid_carbons = total_carbons - 23
    
    # Medium-chain fatty acids have 6-12 carbons
    if fatty_acid_carbons < 6:
        return False, f"Fatty acid chain too short ({fatty_acid_carbons} carbons, need 6-12)"
    if fatty_acid_carbons > 12:
        return False, f"Fatty acid chain too long ({fatty_acid_carbons} carbons, need 6-12)"

    # Check for ribose sugar
    ribose_pattern = Chem.MolFromSmarts("OC1C(O)C(O)C(O1)CN")
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "Missing ribose sugar moiety"

    return True, f"Medium-chain fatty acyl-CoA(4-) with {fatty_acid_carbons} carbons in fatty acid chain"
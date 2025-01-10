"""
Classifies: CHEBI:65111 3-substituted propionyl-CoA(4-)
"""
"""
Classifies: 3-substituted propionyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_substituted_propionyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-substituted propionyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-substituted propionyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for required number of negative charges (4)
    charge_pattern = Chem.MolFromSmarts("[O-]")
    charge_matches = mol.GetSubstructMatches(charge_pattern)
    if len(charge_matches) != 4:
        return False, f"Found {len(charge_matches)} negative charges, need exactly 4"

    # Check for CoA backbone components
    
    # Thioester linkage
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"
    
    # Adenine base
    adenine_pattern = Chem.MolFromSmarts("c1nc(N)c2ncnc2n1")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "No adenine base found"
    
    # Phosphate groups
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)([O-])")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate groups found"
    
    # Pantetheine moiety
    pantetheine_pattern = Chem.MolFromSmarts("NCCC(=O)NCCS")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "No pantetheine moiety found"

    # Check for ribose
    ribose_pattern = Chem.MolFromSmarts("C1C(O)C(O)C(O)C1")
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "No ribose moiety found"

    # Count carbons in acyl chain
    acyl_chain_carbons = 0
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if thioester_matches:
        # Get the carbon atom of the thioester
        thioester_carbon = thioester_matches[0][0]
        # Count carbons in the chain attached to the thioester
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:  # Carbon atom
                if atom.GetIdx() != thioester_carbon:
                    acyl_chain_carbons += 1

    if acyl_chain_carbons < 3:
        return False, "Acyl chain too short for 3-substituted propionyl-CoA"

    return True, "Contains CoA(4-) moiety with appropriate thioester linkage and substituted acyl chain"
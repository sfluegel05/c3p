"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2-enoyl-CoA based on its SMILES string.
    A 2-enoyl-CoA is an unsaturated fatty acyl-CoA in which the S-acyl group 
    contains a double bond between positions 2 and 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-enoyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES with stereochemistry
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA backbone pattern
    # Look for adenine nucleobase
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "No CoA moiety found (missing adenine)"
    
    # Look for thioester group (-C(=O)S-)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester group found"

    # Look for pantetheine part
    pantetheine_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Missing pantetheine part"

    # Look for characteristic alpha-beta unsaturation (position 2-3)
    # Pattern matches both cis and trans configurations
    alpha_beta_pattern = Chem.MolFromSmarts("[C,H]-[CX3]=[CX3]-C(=O)S")
    if not mol.HasSubstructMatch(alpha_beta_pattern):
        return False, "No double bond between positions 2 and 3 relative to thioester"

    # Check for phosphate groups (CoA has 3 phosphates)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)([O,OH])[O,OH]")
    phosphate_matches = len(mol.GetSubstructMatches(phosphate_pattern))
    if phosphate_matches < 3:
        return False, f"Found only {phosphate_matches} phosphate groups, need at least 3"

    # Check for ribose sugar
    ribose_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "Missing or incorrect ribose sugar moiety"

    return True, "Contains CoA moiety and unsaturated fatty acyl group with double bond between positions 2 and 3"
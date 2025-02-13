"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule belongs to the class 11,12-saturated fatty acyl-CoA(4-) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is 11,12-saturated fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to obtain an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern to find a long chain fatty acid, flexible based on saturation and chain length
    # This should match long hydrophobic carbon chains
    fatty_acid_pattern = Chem.MolFromSmarts("C[CX4,CX3,CX2]~[CX4,CX3,CX2]~[CX4,CX3,CX2]~[CX4,CX3,CX2]~[CX4,CX3,CX2]~[CX4,CX3,CX2]~[CX4,CX3,CX2]")
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "No or incomplete long chain fatty acid detected"

    # Recognize the coenzyme A moiety, allowing flexibility with bonds and charges
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OCC1OCC(O)C1OP(=O)([O-])[O-]")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"

    # Inclusion of chiral centers and specific 3-hydroxyl group identification
    hydroxyl_pattern = Chem.MolFromSmarts("[C@@H](O)CC(=O)S")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Chiral hydroxyl group not in expected position"

    return True, "Molecule matches the structural patterns for 11,12-saturated fatty acyl-CoA(4-)"
"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule belongs to the class 11,12-saturated fatty acyl-CoA(4-)
    based on its SMILES string.
    
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

    # More specific pattern for a long chain fatty acid typically seen in 11,12-saturated fatty chains
    fatty_acid_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCCCC(=O)CC(=O)")  # Adjusted generic representation
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "No or incomplete long chain fatty acid detected"

    # Recognize the coenzyme A moiety with necessary functional groups
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OCC(O)C1O[C@H]1")  # Pattern revised for clarity
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"

    # Improved chiral center and specific hydroxyl group identification
    # Check for presence of a 3-hydroxy structure stereocenter
    hydroxyl_pattern = Chem.MolFromSmarts("[C@@H](O)C")  # Simplified and made more relevant
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Chiral hydroxyl group not in expected position"

    return True, "Molecule matches the structural patterns for 11,12-saturated fatty acyl-CoA(4-)"
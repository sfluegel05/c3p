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
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Long chain fatty acid pattern: allow variable lengths, with possible unsaturations
    fatty_acid_pattern = Chem.MolFromSmarts("C(CCCCCCCCCCCCCCC)C")
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "No long chain fatty acid detected"

    # Improved CoA moiety pattern
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])O[C@H]1O[C@H](C)C(O)C1OP(=O)([O-])[O-]")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"

    # Analyze stereochemistry and specific hydroxyl locations
    hydroxyl_pattern = Chem.MolFromSmarts("[C@@H](O)CC(=O)S")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Chiral hydroxyl not in expected position"

    return True, "Molecule matches the structural patterns for 11,12-saturated fatty acyl-CoA(4-)"
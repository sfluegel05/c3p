"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA consists of a 3-hydroxy fatty acid moiety linked to coenzyme A.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    # CoA moiety pattern: Adenine - Ribose - Phosphate - Pantetheine
    coA_pattern = Chem.MolFromSmarts("NC1=NC=NC2=C1N=C(NC2=O)N3[C@H]4O[C@H]([C@@H](O[C@H]4COP(O)(=O)O)O[P@](=O)([O-])O)C=[N+]5(CCSC(=O)C=C)C(=O)C5=O")

    # 3-hydroxy fatty acid pattern: Long carbon chain with an OH on the 3rd carbon
    hydroxy_fatty_acid_pattern = Chem.MolFromSmarts("[CH2][CH2][CH](O)[CH2]")

    # Check for 3-hydroxy fatty acid moiety
    if not mol.HasSubstructMatch(hydroxy_fatty_acid_pattern):
        return False, "No 3-hydroxy fatty acid moiety found"

    # Check for CoA moiety
    if not mol.HasSubstructMatch(coA_pattern):
        return False, "No coenzyme A moiety found"

    return True, "Molecule contains a 3-hydroxy fatty acid moiety linked to coenzyme A, indicative of 3-hydroxy fatty acyl-CoA"
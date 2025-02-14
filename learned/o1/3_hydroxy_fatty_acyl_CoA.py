"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
"""
Classifies: 3-hydroxy fatty acyl-CoA
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA results from the condensation of Coenzyme A (CoA) with a 3-hydroxy fatty acid via a thioester linkage.

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

    # Define SMARTS pattern for pantetheine moiety (part of CoA)
    pantetheine_smarts = 'SCCNC(=O)CCNC(=O)'
    pantetheine_pattern = Chem.MolFromSmarts(pantetheine_smarts)
    if pantetheine_pattern is None:
        return False, "Error parsing pantetheine SMARTS pattern"

    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Pantetheine moiety (part of CoA) not found"

    # Define SMARTS pattern for adenine ring (part of CoA)
    adenine_smarts = 'n1cnc2c(n1)ncnc2'
    adenine_pattern = Chem.MolFromSmarts(adenine_smarts)
    if adenine_pattern is None:
        return False, "Error parsing adenine SMARTS pattern"

    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Adenine moiety (part of CoA) not found"

    # Define SMARTS pattern for 3-hydroxy fatty acyl chain attached via thioester linkage
    # Looking for C(=O)S-C-C(OH)-C, where the OH is on the third carbon from the carbonyl
    thioester_hydroxy_smarts = 'C(=O)SCC[CH](O)'
    thioester_hydroxy_pattern = Chem.MolFromSmarts(thioester_hydroxy_smarts)
    if thioester_hydroxy_pattern is None:
        return False, "Error parsing thioester hydroxy SMARTS pattern"

    if not mol.HasSubstructMatch(thioester_hydroxy_pattern):
        return False, "No 3-hydroxy group on acyl chain attached via thioester linkage found"

    return True, "Contains 3-hydroxy fatty acyl-CoA structure"
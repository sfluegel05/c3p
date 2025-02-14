"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors


def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA has a coenzyme A moiety linked via a thioester to the carbonyl of a 3-hydroxy fatty acid.

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

    # Combined SMARTS pattern for CoA, thioester, and 3-hydroxy fatty acid chain
    # This pattern specifies the CoA, a thioester (-C(=O)S-) a carbon attached to the S, a carbon attached to the previous carbon, and then the 3-hydroxy group attached to this last carbon. Then, any number of carbons should be attached to it
    # The pattern also specifies that there should be at least 3 carbons in the chain
    combined_pattern = Chem.MolFromSmarts("[#6](=[#8])[#16][#6][#6]([#8])[#6]~[#6]~[#6]~[#6]COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)n1cnc2c(N)ncnc12")
  
    if mol.HasSubstructMatch(combined_pattern):
        return True, "Contains CoA, thioester, and a 3-hydroxy fatty acid chain."
    else:
        return False, "Molecule does not match the required pattern for 3-hydroxy fatty acyl-CoA"
"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2-enoyl-CoA based on its SMILES string.
    A 2-enoyl-CoA has an unsaturated fatty acyl-CoA with a double bond between positions 2 and 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-enoyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for 2-enoyl characteristic C=C bond between positions 2 and 3
    enoyl_pattern = Chem.MolFromSmarts("C/C=C/C(=O)")
    if not mol.HasSubstructMatch(enoyl_pattern):
        return False, "No 2-enoyl C=C bond found between positions 2 and 3"
    
    # Check for CoA linkage pattern
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA linkage pattern not found"

    return True, "Contains 2-enoyl C=C bond with CoA linkage"
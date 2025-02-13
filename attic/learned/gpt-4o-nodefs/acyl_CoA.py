"""
Classifies: CHEBI:17984 acyl-CoA
"""
from rdkit import Chem

def is_acyl_CoA(smiles: str) -> (bool, str):
    """
    Classifies the molecule as acyl-CoA based on SMILES string.
    Acyl-CoA is characterized by a fatty acid chain bound to Coenzyme A via a thioester bond.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as acyl-CoA, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern to identify CoA portion
    coa_pattern = Chem.MolFromSmarts("C(=O)NCCSC(=O)C[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing CoA moiety"

    # Pattern to identify thioester bond, a key feature of acyl-CoA
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCCN")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester bond with CoA moiety"

    return True, "Contains CoA moiety with a thioester bond, characteristic of acyl-CoA"
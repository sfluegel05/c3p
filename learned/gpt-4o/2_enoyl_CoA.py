"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
from rdkit import Chem

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2-enoyl-CoA based on its SMILES string.
    A 2-enoyl-CoA is characterized by an unsaturated fatty acyl-CoA 
    with a double bond between carbon positions 2 and 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-enoyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES to Molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Finding the 2-enoyl pattern involving C2=C3 with carbon chain pattern
    # This pattern captures the double bond between positions 2 and 3 adjacent to the carbonyl
    enoyl_pattern = Chem.MolFromSmarts("CC=C(=O)C")
    if not mol.HasSubstructMatch(enoyl_pattern):
        return False, "No characteristic C2=C3 enoyl C=C bond found"

    # Standard CoA linkage check
    # Broader and more flexible pattern for CoA linkage
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H](O)C(O)C(O)C1OP1")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA linkage pattern not found"

    return True, "Contains 2-enoyl C2=C3 bond with CoA linkage"
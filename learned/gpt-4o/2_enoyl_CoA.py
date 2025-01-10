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

    # Finding the enoyl pattern involving C2=C3
    # Using an adjusted pattern to match more variants of 2-enoyl including stereochemistry
    enoyl_pattern = Chem.MolFromSmarts("CC=CC(=O)S")  # Simplified for core detection
    if not mol.HasSubstructMatch(enoyl_pattern):
        return False, "No characteristic C2=C3 enoyl C=C bond found"

    # Standard CoA linkage check
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(O)(=O)OP(O)(=O)OC1C(O)C(O)C1O")  # Generalized a bit
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA linkage pattern not found"

    return True, "Contains 2-enoyl C2=C3 bond with CoA linkage"
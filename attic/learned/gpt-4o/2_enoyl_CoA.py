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

    # Refined searching for 2-enoyl pattern
    # This updated pattern attempts to better match the exemplified enoyl groups
    enoyl_pattern = Chem.MolFromSmarts("[#6]=[#6][CX3](=O)")
    if not mol.HasSubstructMatch(enoyl_pattern):
        return False, "No characteristic enoyl C=C bond found near carbonyl group"

    # Refined Standard CoA linkage checking for more general CoA patterns
    # Adding flexibility to handle different representations of CoA linkages
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[CX4][OX2]")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA linkage pattern not found"

    return True, "Contains 2-enoyl C=C bond with CoA linkage"
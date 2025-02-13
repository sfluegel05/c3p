"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
from rdkit import Chem

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acyl-CoA based on its SMILES string.
    A fatty acyl-CoA results from the condensation of the thiol group of coenzyme A with
    the carboxy group of an unsaturated fatty acid.

    Args:
    - smiles (str): SMILES string of the molecule

    Returns:
    - bool: True if molecule is an unsaturated fatty acyl-CoA, False otherwise
    - str: Reason for classification
    """
    # Parse the SMILES string into an RDKit mol object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the thioester linkage: C(=O)S
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Check for the Coenzyme A portion using a representative pattern
    coenzymeA_pattern = Chem.MolFromSmarts('C(=O)NC(C)(C)COP(O)(O)=O')
    if not mol.HasSubstructMatch(coenzymeA_pattern):
        return False, "Coenzyme A moiety pattern not found"

    # Check for unsaturated fatty acid chain (presence of C=C)
    unsaturated_fatty_acid_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(unsaturated_fatty_acid_pattern):
        return False, "No unsaturation (C=C bond) found"

    return True, "Contains both unsaturated fatty acid and CoA moiety"
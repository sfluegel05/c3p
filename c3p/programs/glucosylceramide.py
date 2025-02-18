"""
Classifies: CHEBI:36500 glucosylceramide
"""
from rdkit import Chem

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    A glucosylceramide contains a glucose moiety, a sphingosine backbone, and a long fatty acid chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosylceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify β-D-glucose moiety with possible stereochemistry
    glucose_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H](O)[C@H](CO)O[C@H]1")
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No -glucose group found"

    # Identify amide linkage (N-C(=O))
    amide_pattern = Chem.MolFromSmarts("NC(=O)")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage found"

    # Identify sphingosine backbone (typical C18 chain with amine and alcohol functionalities)
    sphingosine_pattern = Chem.MolFromSmarts("[#6]-[#6]-[#6]-[#6]-[#6](=[#8])-[#7]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](-[#6])-[#6]")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone found"

    # Identify long fatty acid chain (carbon chain extending from backbone)
    long_chain_pattern = Chem.MolFromSmarts("C" * 16)  # At least 16 carbons
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long fatty acid chain found"

    return True, "Contains glucose moiety linked to sphingosine backbone with a fatty acid chain"
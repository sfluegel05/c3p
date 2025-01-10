"""
Classifies: CHEBI:17297 UDP-sugar
"""
from rdkit import Chem

def is_UDP_sugar(smiles: str):
    """
    Determine if a molecule is a UDP-sugar based on its SMILES string. A UDP-sugar is a pyrimidine nucleotide-sugar with UDP as the nucleotide component attached via an anomeric diphosphate linkage.

    Args:
        smiles (str): The SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a UDP-sugar, False otherwise.
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Redefine uridine pattern with broader flexibility
    uridine_pattern = Chem.MolFromSmarts("C1=CN(C(=O)N=C1)C2C(C(C(O2)CO)O)O")
    # Checking uridine presence
    if not mol.HasSubstructMatch(uridine_pattern):
        return False, "Uridine moiety not found"

    # Redefine a more general diphosphate pattern
    diphosphate_pattern = Chem.MolFromSmarts("OP(O)(=O)OP(O)(=O)O")
    # Checking diphosphate linkage
    if not mol.HasSubstructMatch(diphosphate_pattern):
        return False, "Diphosphate linkage not found"

    return True, "Contains UDP moiety with a sugar via diphosphate linkage"
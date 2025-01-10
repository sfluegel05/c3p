"""
Classifies: CHEBI:27325 xanthophyll
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    Xanthophylls are characterized by long conjugated carbon chains and the presence of oxygen functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a xanthophyll, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for long conjugated carbon chain
    conjugated_chain_pattern = Chem.MolFromSmarts("[C]=[C]-[C]=[C]-[C]=[C]")  # Pattern for a chain with alternating single and double bonds
    if not mol.HasSubstructMatch(conjugated_chain_pattern):
        return False, "No long conjugated carbon chain found, typical in xanthophylls"
        
    # Check for presence of oxygen functionalities (-OH, =O, epoxy)
    contains_oxygen = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # Atomic number 8 is Oxygen
            contains_oxygen = True
            break
    if not contains_oxygen:
        return False, "No oxygen functionality present, needed for classification as a xanthophyll"

    return True, "Contains long conjugated carbon chain and oxygen functionality typical of xanthophylls"
"""
Classifies: CHEBI:36092 clavulone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a clavulone, False otherwise
        str: Reason for classification
    """
     
    # Parse SMILES into RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for aliphatic chain pattern
    aliphatic_chain = Chem.MolFromSmarts("C/C=C\\CCCCC")
    if not mol.HasSubstructMatch(aliphatic_chain):
        return False, "Missing long aliphatic chain"

    # Check for presence of halogens
    halogen = Chem.MolFromSmarts("[Cl,Br,I]")
    if mol.HasSubstructMatch(halogen):
        halogen_present = True
    else:
        halogen_present = False

    # Check for ester groups
    ester_pattern = Chem.MolFromSmarts("OC(=O)C")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "Missing ester group"

    # All characteristic patterns found
    if halogen_present:
        return True, "Contains long aliphatic chain, ester group, and halogen"
    else:
        return True, "Contains long aliphatic chain and ester group"

    return None, None
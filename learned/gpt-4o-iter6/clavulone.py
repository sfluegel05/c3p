"""
Classifies: CHEBI:36092 clavulone
"""
from rdkit import Chem

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone based on its SMILES string.
    A clavulone is derived from marine corals and is classified as an esterified prostanoid.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a clavulone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for ester group pattern
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group found"
    
    # Look for conjugated double bond pattern in a ring
    conjugated_diene_pattern = Chem.MolFromSmarts("C-C=C-C=C") 
    if not mol.HasSubstructMatch(conjugated_diene_pattern):
        return False, "No conjugated diene found"
    
    # Check for presence of halogen
    halogen_pattern = Chem.MolFromSmarts("[Cl,Br,I]")
    if not mol.HasSubstructMatch(halogen_pattern):
        return False, "No halogen found"

    # Additional checks for long carbon chain (from marine origin)
    long_chain_pattern = Chem.MolFromSmarts("C=C-C=C-C=C-C=C")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long carbon chain found"
    
    return True, "Matches clavulone structure"

# Example call
# result, reason = is_clavulone("ClC=1C(=O)[C@@]([C@@](O)(C/C=C\\CCCCC)C1)([C@@H](OC(=O)C)[C@H](OC(=O)C)[C@@H](OC(=O)C)CCCC(OC)=O)[H]")
# print(result, reason)
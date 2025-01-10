"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2-enoyl-CoA derivative based on its SMILES string.

    A 2-enoyl-CoA typically features a CoA moiety connected by a thioester linkage to an acyl chain
    with an enoyl group (C=C) towards the beginning of the chain.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-enoyl-CoA derivative, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # CoA moiety pattern (a simplified pattern can be used)
    coa_pattern = Chem.MolFromSmarts("C1=NC=NC2=NC=CN=C21")  # Simplified pattern for purine-like structure 
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"
    
    # Enoyl pattern (with a double bond near start of chain)
    enoyl_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(enoyl_pattern):
        return False, "No enoyl group found"
    
    # Thioester linkage pattern
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    return True, "Contains CoA moiety, enoyl group, and thioester linkage"
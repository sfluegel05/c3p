"""
Classifies: CHEBI:17984 acyl-CoA
"""
"""
Classifies: CHEBI:17996 acyl-CoA
"""
from rdkit import Chem

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    An acyl-CoA is a thioester formed between coenzyme A's thiol group and a carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the key substructure: thioester (S-C(=O)) connected to CoA's pantetheine chain
    # Pattern matches S connected to both the acyl group (C=O) and the CoA backbone (CCNC...)
    acyl_coa_pattern = Chem.MolFromSmarts("[S](C(=O))CCNC(=O)CCNC(=O)")
    
    if mol.HasSubstructMatch(acyl_coa_pattern):
        return True, "Contains thioester group linked to CoA pantetheine chain"
    else:
        return False, "Missing thioester linkage to CoA backbone"
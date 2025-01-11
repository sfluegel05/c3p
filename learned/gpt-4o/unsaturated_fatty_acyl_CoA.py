"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
from rdkit import Chem

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acyl-CoA based on its SMILES string.
    An unsaturated fatty acyl-CoA contains features of an unsaturated fatty acid bonded to coenzyme A.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of unsaturated fatty acid part (C=C double bonds)
    unsaturated_fatty_acid_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(unsaturated_fatty_acid_pattern):
        return False, "No unsaturated carbon-carbon double bonds found"

    # Check for the presence of coenzyme A substructure
    # Coenzyme A contains a pantetheine tail and terminal phosphate groups; use key markers like C(=O)S (thioester),
    # and phosphate group patterns (OP(=O)(O)O)
    coa_thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    coa_phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O")
    
    if not mol.HasSubstructMatch(coa_thioester_pattern):
        return False, "Thioester linkage not found, indicating absence of CoA attachment"

    if not mol.HasSubstructMatch(coa_phosphate_pattern):
        return False, "Phosphate groups typical of CoA not found"

    return True, "Contains unsaturated fatty acid chain and coenzyme A substructure with thioester linkage"
"""
Classifies: CHEBI:17984 acyl-CoA
"""
"""
Classifies: CHEBI:37577 acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    An acyl-CoA is a thioester formed from the condensation of the thiol group of coenzyme A
    with the carboxy group of any carboxylic acid.

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

    # Define SMARTS pattern for the Coenzyme A moiety (simplified key features)
    coenzymeA_smarts = "NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(O)[C@@H](OP(=O)(O)[O-])[C@H](O)[C@H]1O"
    coenzymeA_pattern = Chem.MolFromSmarts(coenzymeA_smarts)
    if coenzymeA_pattern is None:
        return False, "Failed to parse Coenzyme A SMARTS pattern"

    # Check for Coenzyme A moiety
    if not mol.HasSubstructMatch(coenzymeA_pattern):
        return False, "Coenzyme A moiety not found"

    # Define SMARTS pattern for the thioester linkage (S-C(=O)-C)
    thioester_smarts = "SC(=O)C"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    if thioester_pattern is None:
        return False, "Failed to parse thioester SMARTS pattern"

    # Check for thioester linkage
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester linkage not found"

    # Ensure the acyl group is attached via thioester bond
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester linkage not found"
    else:
        # Verify that the sulfur is connected to the CoA moiety
        sulfur_idx = thioester_matches[0][0]
        # Check if sulfur is part of the Coenzyme A moiety
        coa_match = mol.GetSubstructMatch(coenzymeA_pattern)
        if sulfur_idx not in coa_match:
            return False, "Thioester linkage is not connected to Coenzyme A"

    return True, "Contains Coenzyme A moiety linked via thioester bond to an acyl group"
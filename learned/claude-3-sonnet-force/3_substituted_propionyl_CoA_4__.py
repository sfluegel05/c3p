"""
Classifies: CHEBI:65111 3-substituted propionyl-CoA(4-)
"""
"""
Classifies: CHEBI:84547 3-substituted propionyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_substituted_propionyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-substituted propionyl-CoA(4-) based on its SMILES string.
    A 3-substituted propionyl-CoA(4-) is an acyl-CoA(4-) oxoanion arising from deprotonation
    of the phosphate and diphosphate OH groups of any 3-substituted propionyl-CoA.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-substituted propionyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for CoA backbone
    coa_backbone_pattern = Chem.MolFromSmarts("C(C(C(=O)NCCC(=O)NCCSC)O)O[P@@](=O)(O)[O-].[n&r5]1[c&r6]2[c&r7]([n&r5]1)[n&r6][c&r7][n&r6][c&r7]2[N&r6]")
    if not mol.HasSubstructMatch(coa_backbone_pattern):
        return False, "Missing CoA backbone"
    
    # Look for 3-substituted propionyl group attached to sulfur
    propionyl_pattern = Chem.MolFromSmarts("S(=O)(=O)CC(C)C")
    propionyl_matches = mol.GetSubstructMatches(propionyl_pattern)
    if not propionyl_matches:
        return False, "Missing 3-substituted propionyl group attached to sulfur"
    
    # Check for double bonds in the acyl chain
    acyl_chain_pattern = Chem.MolFromSmarts("S(=O)(=O)CCC=C")
    acyl_chain_matches = mol.GetSubstructMatches(acyl_chain_pattern)
    
    # Look for deprotonated phosphate and diphosphate groups
    phosphate_pattern = Chem.MolFromSmarts("[O-,P]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 3:
        return False, "Missing deprotonated phosphate and diphosphate groups"
    
    if acyl_chain_matches:
        return True, "Contains CoA backbone with 3-substituted propionyl group, double bonds in acyl chain, and deprotonated phosphate/diphosphate groups"
    else:
        return True, "Contains CoA backbone with 3-substituted propionyl group and deprotonated phosphate/diphosphate groups"
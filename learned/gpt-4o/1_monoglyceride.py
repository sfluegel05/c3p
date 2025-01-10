"""
Classifies: CHEBI:35759 1-monoglyceride
"""
from rdkit import Chem

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride is a glycerol moiety with a single acyl group attached at
    the primary position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 1-monoglyceride, False otherwise
        str: Explanation of the result
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone (2 primary alcohols and 1 secondary alcohol)
    glycerol_backbone = Chem.MolFromSmarts("C(O)C(O)CO")
    if not mol.HasSubstructMatch(glycerol_backbone):
        return False, "No glycerol backbone found"

    # Check for ester group at primary position (-O-C(=O)-)
    ester_primary_position = Chem.MolFromSmarts("[CH2]OC(=O)[*]")
    if not mol.HasSubstructMatch(ester_primary_position):
        return False, "No ester linkage at primary position found"

    # Confirm that there is only one ester
    ester_group_check = mol.GetSubstructMatches(Chem.MolFromSmarts("COC(=O)[*]"))
    if len(ester_group_check) != 1:
        return False, f"Found {len(ester_group_check)} ester groups, need exactly 1"

    # Consider chiral center criteria (i.e., sn-glycerol)
    chiral_sn_glycerol = Chem.MolFromSmarts("[C@@H]([O])[C@H](O)")
    if not mol.HasSubstructMatch(chiral_sn_glycerol):
        return False, "Specific stereochemistry for sn-glycerol not found"

    return True, "Contains glycerol backbone with a single acyl group esterified at a primary position"
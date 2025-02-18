"""
Classifies: CHEBI:17962 CDP-diacylglycerol
"""
"""
Classifies: CHEBI:17855 CDP-diacylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_CDP_diacylglycerol(smiles: str):
    """
    Determines if a molecule is a CDP-diacylglycerol based on its SMILES string.
    A CDP-diacylglycerol consists of a glycerol backbone with two acyl groups at positions 1 and 2,
    and a cytidine diphosphate (CDP) group at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a CDP-diacylglycerol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly two ester groups (acyl chains at positions 1 and 2)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Verify glycerol backbone (C1-C2-C3 where C1 and C2 have ester groups)
    # Modified to account for stereochemistry and branching
    glycerol_backbone = Chem.MolFromSmarts("[CH2]-[CH](-O-C(=O)-*)-[CH2]")
    if not mol.HasSubstructMatch(glycerol_backbone):
        return False, "Glycerol backbone with ester groups not found"

    # Check diphosphate bridge (P-O-P-O connection)
    # More flexible pattern to match different protonation states
    diphosphate_pattern = Chem.MolFromSmarts("P-O-P")
    if not mol.HasSubstructMatch(diphosphate_pattern):
        return False, "Diphosphate bridge not found"

    # Updated cytosine pattern to match actual examples (n1c(=O)nc(N)cc1)
    cytosine_pattern = Chem.MolFromSmarts("n1c(=O)nc(N)cc1")
    if not mol.HasSubstructMatch(cytosine_pattern):
        return False, "Cytosine base not detected"

    # Simplified ribose check (five-membered ring with O and hydroxyls)
    ribose_pattern = Chem.MolFromSmarts("[C@H]1O[C@H](CO)[C@H](O)[C@@H]1O")
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "Ribose sugar not detected"

    # Verify connectivity between components:
    # Glycerol -> diphosphate -> ribose -> cytosine
    # Using a more robust connectivity check
    cdp_pattern = Chem.MolFromSmarts(
        "[CH2]-[CH](-O-C(=O)-*)-[CH2]-O-P-O-P-O-C1-C-O-[C@H]2O[C@H](CO)[C@H](O)[C@@H]2O.n3c(=O)nc(N)cc3"
    )
    if not mol.HasSubstructMatch(cdp_pattern):
        return False, "CDP-glycerol connectivity not confirmed"

    return True, "Contains glycerol with two acyl groups and properly connected CDP group"
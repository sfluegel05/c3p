"""
Classifies: CHEBI:17297 UDP-sugar
"""
"""
Classifies: CHEBI:15839 UDP-sugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_UDP_sugar(smiles: str):
    """
    Determines if a molecule is a UDP-sugar based on its SMILES string.
    A UDP-sugar is a pyrimidine nucleotide-sugar having UDP as the nucleotide component attached to an unspecified sugar via an anomeric diphosphate linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a UDP-sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the uracil ring pattern (pyrimidine base)
    uracil_pattern = Chem.MolFromSmarts("[nH0,nH1]1ccc(=O)[nH0,nH1]c1=O")
    if not mol.HasSubstructMatch(uracil_pattern):
        return False, "No uracil ring found"

    # Look for the ribose sugar attached to uracil - more general pattern
    ribose_pattern = Chem.MolFromSmarts("[nH0,nH1]1ccc(=O)[nH0,nH1]c1=O[C@H]1[C@H](O)[C@H](O)[C@H](O)CO1")
    if not mol.HasSubstructMatch(ribose_pattern):
        # Relax the pattern to allow for variations in the ribose structure
        ribose_pattern = Chem.MolFromSmarts("[nH0,nH1]1ccc(=O)[nH0,nH1]c1=O[C@H]1[C@H](O)[C@H](O)[C@H](O)C1")
        if not mol.HasSubstructMatch(ribose_pattern):
            return False, "No ribose sugar attached to uracil found"

    # Look for the diphosphate group (two phosphate groups in sequence)
    diphosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])OP(=O)([OX2])[OX2]")
    diphosphate_matches = mol.GetSubstructMatches(diphosphate_pattern)
    if len(diphosphate_matches) == 0:
        return False, "No diphosphate group found"

    # Check if the diphosphate is linked to a sugar (any sugar)
    # We assume that the sugar is any ring structure with multiple hydroxyl groups
    sugar_pattern = Chem.MolFromSmarts("[C@H]1([C@H](O)[C@H](O)[C@H](O)[C@H](O1))")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) == 0:
        return False, "No sugar moiety found"

    # Check if the diphosphate is linked to the sugar
    diphosphate_atoms = set()
    for match in diphosphate_matches:
        diphosphate_atoms.update(match)
    
    sugar_atoms = set()
    for match in sugar_matches:
        sugar_atoms.update(match)
    
    if not diphosphate_atoms.intersection(sugar_atoms):
        return False, "Diphosphate not linked to sugar"

    return True, "Contains UDP moiety linked to a sugar via an anomeric diphosphate linkage"
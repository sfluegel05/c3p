"""
Classifies: CHEBI:24128 furanocoumarin
"""
"""
Classifies: CHEBI:32798 furanocoumarin

A furanocoumarin is any furochromene that consists of a furan ring fused with a coumarin.
The fusion may occur in different ways to give several isomers.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a furan and a coumarin ring
    furan_pattern = Chem.MolFromSmarts("[o,O]1cccc1")
    coumarin_pattern = Chem.MolFromSmarts("c1ccc2c(c1)cc(=O)oc2")
    if not mol.HasSubstructMatch(furan_pattern) or not mol.HasSubstructMatch(coumarin_pattern):
        return False, "Missing furan or coumarin ring"

    # Check for the presence of a fused ring system
    fused_ring_count = rdMolDescriptors.CalcNumAromaticRings(mol) + rdMolDescriptors.CalcNumAliphaticRings(mol) - mol.GetRingInfo().NumRings()
    if fused_ring_count < 1:
        return False, "No fused ring system found"

    # Check for the presence of ether or ester groups (common in furanocoumarins)
    ether_pattern = Chem.MolFromSmarts("[OX2]")
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    if not (mol.HasSubstructMatch(ether_pattern) or mol.HasSubstructMatch(ester_pattern)):
        return False, "No ether or ester groups found"

    # Check molecular weight range (typical for furanocoumarins)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 600:
        return False, "Molecular weight outside typical range for furanocoumarins"

    return True, "Contains furanocoumarin scaffold (fused furan and coumarin rings, ether/ester groups, reasonable molecular weight)"
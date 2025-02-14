"""
Classifies: CHEBI:24128 furanocoumarin
"""
"""
Classifies: CHEBI:27727 furanocoumarin

A furanocoumarin is defined as any furochromene that consists of a furan ring fused with a coumarin.
The fusion may occur in different ways to give several isomers.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_furanocoumarin(smiles: str) -> tuple[bool, str]:
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

    # Check for the presence of a coumarin ring system
    coumarin_pattern = Chem.MolFromSmarts("c1cc2ccoc2cc1=O")
    if not mol.HasSubstructMatch(coumarin_pattern):
        return False, "No coumarin ring system found"

    # Check for the presence of a furan ring
    furan_pattern = Chem.MolFromSmarts("c1occc1")
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring found"

    # Check if the coumarin and furan rings are fused
    fused_system_pattern = Chem.MolFromSmarts("c1cc2ccoc2cc1oc3cccc3")
    if not mol.HasSubstructMatch(fused_system_pattern):
        return False, "Coumarin and furan rings are not fused"

    # Additional checks for furanocoumarin characteristics
    ring_info = mol.GetRingInfo()
    fused_ring_count = len([r for r in ring_info.AtomRings() if len(r) > 6])
    if fused_ring_count != 1:
        return False, "Incorrect number of fused rings"

    # Check for appropriate molecular weight range
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 200 or mol_weight > 500:
        return False, "Molecular weight outside the typical range for furanocoumarins"

    # Check for the presence of oxygen atoms (typically 3-4)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 3 or oxygen_count > 4:
        return False, "Unexpected number of oxygen atoms for a furanocoumarin"

    return True, "Contains a fused coumarin and furan ring system (furanocoumarin)"
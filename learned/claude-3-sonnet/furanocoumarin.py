"""
Classifies: CHEBI:24128 furanocoumarin
"""
"""
Classifies: CHEBI:51287 furanocoumarin

A furanocoumarin is any furochromene that consists of a furan ring fused with a coumarin.
The fusion may occur in different ways to give several isomers.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
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

    # Check for furan ring
    furan_pattern = Chem.MolFromSmarts("O1C=CC=C1")
    furan_match = mol.HasSubstructMatch(furan_pattern)
    if not furan_match:
        return False, "No furan ring found"

    # Check for coumarin ring
    coumarin_pattern = Chem.MolFromSmarts("O=C1C=CC=CC2=C1OCC=C2")
    coumarin_match = mol.HasSubstructMatch(coumarin_pattern)
    if not coumarin_match:
        return False, "No coumarin ring found"

    # Check for fused ring system
    ring_info = mol.GetRingInfo()
    fused_rings = any(len(x) > 1 for x in ring_info.BondRings())
    if not fused_rings:
        return False, "No fused ring system found"

    # Check for ether/ester groups (optional)
    ether_pattern = Chem.MolFromSmarts("[OX2]")
    ether_match = mol.HasSubstructMatch(ether_pattern)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_match = mol.HasSubstructMatch(ester_pattern)

    # Check molecular weight range (optional)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 600:
        return False, "Molecular weight outside typical range for furanocoumarins"

    if ether_match or ester_match:
        return True, "Contains furan and coumarin rings in a fused system, with ether/ester groups"
    else:
        return True, "Contains furan and coumarin rings in a fused system"
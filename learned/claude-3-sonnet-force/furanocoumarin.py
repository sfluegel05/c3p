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

    # Identify the coumarin core structure
    coumarin_pattern = Chem.MolFromSmarts("[O=2]=[c,C]1[c,C][c,C][c,C][c,C][c,C]1")
    coumarin_matches = mol.GetSubstructMatches(coumarin_pattern)

    # Identify the furan ring
    furan_pattern = Chem.MolFromSmarts("c1occc1")
    furan_matches = mol.GetSubstructMatches(furan_pattern)

    # Check if the coumarin and furan rings are fused
    fused_systems = []
    for coumarin_match in coumarin_matches:
        for furan_match in furan_matches:
            common_atoms = set(coumarin_match) & set(furan_match)
            if common_atoms:
                fused_system = list(set(coumarin_match) | set(furan_match))
                fused_systems.append(fused_system)

    if not fused_systems:
        return False, "Coumarin and furan rings are not fused"

    # Additional checks for furanocoumarin characteristics
    for fused_system in fused_systems:
        ring_info = mol.GetRingInfo()
        fused_ring = ring_info.IsFusedRing(fused_system)
        if not fused_ring:
            continue

        # Check for appropriate atom connectivity and bond orders
        coumarin_atoms = [mol.GetAtomWithIdx(idx) for idx in fused_system if idx in coumarin_match]
        furan_atoms = [mol.GetAtomWithIdx(idx) for idx in fused_system if idx in furan_match]

        coumarin_valid = all(atom.GetAtomicNum() in [6, 8] and atom.GetTotalDegree() in [2, 3] for atom in coumarin_atoms)
        furan_valid = all(atom.GetAtomicNum() in [6, 8] and atom.GetTotalDegree() in [2, 3] for atom in furan_atoms)

        if coumarin_valid and furan_valid:
            return True, "Contains a fused coumarin and furan ring system (furanocoumarin)"

    return False, "Does not match the structural characteristics of a furanocoumarin"
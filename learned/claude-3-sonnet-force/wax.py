"""
Classifies: CHEBI:73702 wax
"""
"""
Classifies: CHEBI:26665 wax
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from typing import Tuple

def is_wax(smiles: str) -> Tuple[bool, str]:
    """
    Determines if a molecule is a wax based on its SMILES string.
    Waxes are defined as organic compounds or mixtures of compounds composed of long-chain
    molecules and malleable at ambient temperatures.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for long carbon chains (>= 20 carbons)
    long_chains = False
    for chain in Chem.Mol.GetRingInfo(mol).AtomRings():
        if sum(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in chain) >= 20:
            long_chains = True
            break

    if not long_chains:
        return False, "No long carbon chains (>= 20 carbons) found"

    # Check for common wax functional groups: esters, ethers, alcohols
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ether_pattern = Chem.MolFromSmarts("[OX2][CX4]")
    alcohol_pattern = Chem.MolFromSmarts("[OX1H]")

    has_wax_groups = bool(mol.HasSubstructMatch(ester_pattern) or
                          mol.HasSubstructMatch(ether_pattern) or
                          mol.HasSubstructMatch(alcohol_pattern))

    if not has_wax_groups:
        return False, "No characteristic wax functional groups (esters, ethers, alcohols) found"

    # Check molecular weight - waxes typically > 500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for wax (< 500 Da)"

    # Count rotatable bonds - waxes typically have > 15
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 15:
        return False, "Too few rotatable bonds for wax (< 15)"

    return True, "Contains long carbon chains and characteristic wax functional groups"
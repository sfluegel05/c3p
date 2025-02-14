"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
"""
Classifies: CHEBI:35468 fatty acid methyl ester
"""
from typing import Tuple
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_fatty_acid_methyl_ester(smiles: str) -> Tuple[bool, str]:
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    A fatty acid methyl ester is a fatty acid ester obtained by condensation with methanol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid methyl ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for methyl ester group (-C(=O)-O-C)
    methyl_ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    if not mol.HasSubstructMatch(methyl_ester_pattern):
        return False, "No methyl ester group found"

    # Look for long aliphatic chain
    aliphatic_chain_pattern = Chem.MolFromSmarts("[CH3,CH2]~[CH2]~[CH2,CH3]")
    aliphatic_chain_matches = mol.GetSubstructMatches(aliphatic_chain_pattern)
    if not aliphatic_chain_matches:
        return False, "No aliphatic chain found"

    # Check for uncommon functional groups or heteroatoms
    uncommon_groups = "[!C!H!O]"  # Exclude C, H, O
    uncommon_groups_pattern = Chem.MolFromSmarts(uncommon_groups)
    if mol.HasSubstructMatch(uncommon_groups_pattern):
        return False, "Uncommon functional groups or heteroatoms found"

    # Check for absence of aromatic rings
    aromatic_atoms = [atom.GetIsAromatic() for atom in mol.GetAtoms()]
    if any(aromatic_atoms):
        return False, "Aromatic rings found"

    # Count rotatable bonds to verify aliphatic chain length
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Aliphatic chain too short"

    return True, "Contains a methyl ester group and a long aliphatic chain, without uncommon functional groups or aromatic rings"
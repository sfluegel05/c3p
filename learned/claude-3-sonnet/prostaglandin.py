"""
Classifies: CHEBI:26333 prostaglandin
"""
"""
Classifies: CHEBI:35493 prostaglandin

Prostaglandins are naturally occurring compounds derived from the parent C20 acid, prostanoic acid.
They share a core structure consisting of a cyclopentane ring fused or joined to a cyclopentene ring,
with various substituents, functional groups, and side chains.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for core cyclopentane and cyclopentene rings
    core_pattern = Chem.MolFromSmarts("[C&R1&r5][C&R1&r5][C&R1&r5][C&R1&r5][C&R1&r5]1[C&R2&r5][C&R2&r5]2[C&R2&r5][C&R2&r5][C&R2&r5]12")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Missing core cyclopentane and cyclopentene rings"

    # Look for common functional groups
    hydroxy_pattern = Chem.MolFromSmarts("[OH]")
    keto_pattern = Chem.MolFromSmarts("[C=O]")
    carboxyl_pattern = Chem.MolFromSmarts("[C$(C=O)O]")
    ether_pattern = Chem.MolFromSmarts("[O]")  # Include ethers
    ester_pattern = Chem.MolFromSmarts("[C$(C=O)OC]")  # Include esters
    epoxide_pattern = Chem.MolFromSmarts("[C1OC1]")  # Include epoxides

    functional_groups = [hydroxy_pattern, keto_pattern, carboxyl_pattern, ether_pattern, ester_pattern, epoxide_pattern]
    functional_group_count = sum(len(mol.GetSubstructMatches(pattern)) for pattern in functional_groups)
    if functional_group_count < 2:
        return False, "Insufficient functional groups"

    # Look for lipid chain (long carbon chain)
    lipid_chain_pattern = Chem.MolFromSmarts("[C&H3]~[C&H2]~[C&H2]~[C&H2]~[C&H2]~[C&H2]~[C&H2]")
    if not mol.HasSubstructMatch(lipid_chain_pattern):
        return False, "Missing lipid chain"

    # Check molecular weight and atom counts
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if mol_wt < 300 or mol_wt > 600:
        return False, "Molecular weight out of typical range for prostaglandins"
    if c_count < 15 or c_count > 30:
        return False, "Carbon count out of typical range for prostaglandins"
    if o_count < 2 or o_count > 8:
        return False, "Oxygen count out of typical range for prostaglandins"

    return True, "Molecule contains core cyclopentane and cyclopentene rings, common functional groups, and a lipid chain"
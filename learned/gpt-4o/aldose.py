"""
Classifies: CHEBI:15693 aldose
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aldose(smiles: str) -> (bool, str):
    """
    Determines if a molecule is an aldose based on its SMILES string.
    Aldoses are polyhydroxy aldehydes (H[CH(OH)]nC(=O)H, n >= 2) and their intramolecular hemiacetals.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldose, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check to determine if the structure is an open-chain aldose:
    # Look for an aldehyde group -C(=O)H
    aldehyde_pattern = Chem.MolFromSmarts("[C;H1](=O)[CH2,CH0]")
    if mol.HasSubstructMatch(aldehyde_pattern):
        # Expected polyhydroxy structure (-OH groups)
        hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
        hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
        if len(hydroxyl_matches) < 2:
            return False, f"Not enough hydroxyl groups found, got {len(hydroxyl_matches)}"

        carbon_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
        if carbon_count < 3:
            return False, f"Insufficient carbon backbone, found {carbon_count}"

        return True, "Structure is consistent with open-chain form of an aldose"

    # Determine if it is a cyclic form with hemiacetal/ketal form:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        # Look for ether linkage with adjacent hydroxyl (O-C-O-H) which indicates hemiacetal
        ether_pattern = Chem.MolFromSmarts("[C](O)([CH2,CH0])O")
        if mol.HasSubstructMatch(ether_pattern):
            # Check for adequate carbon and hydroxyl groups in the cyclic form
            carbon_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
            hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
            if carbon_count >= 3 and len(hydroxyl_matches) >= 2:
                return True, "Structure is consistent with cyclic hemiacetal form of an aldose"

    return False, "Structure does not fit criteria for an aldose"
"""
Classifies: CHEBI:64611 ether lipid
"""
"""
Classifies: Ether lipid (CHEBI: ???)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    Ether lipids have a glycerol backbone with at least one ether-linked alkyl chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ether lipid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for glycerol backbone (C-C-C with three oxygen connections)
    glycerol = Chem.MolFromSmarts("[CH2]-[CH]-[CH2]")
    if not mol.HasSubstructMatch(glycerol):
        return False, "No glycerol backbone"

    # Find oxygen atoms connected to two carbons (ether linkages)
    ether_pattern = Chem.MolFromSmarts("[C][O][C]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    if not ether_matches:
        return False, "No ether linkages found"

    # Check if any ether oxygen is attached to the glycerol
    glycerol_atoms = set(mol.GetSubstructMatch(glycerol))
    ether_attached = False
    for match in ether_matches:
        oxygen_idx = match[1]
        oxygen = mol.GetAtomWithIdx(oxygen_idx)
        neighbors = [n.GetIdx() for n in oxygen.GetNeighbors()]
        # Check if one neighbor is in glycerol and the other is a carbon chain
        if len(neighbors) == 2:
            if (neighbors[0] in glycerol_atoms) or (neighbors[1] in glycerol_atoms):
                ether_attached = True
                break

    if not ether_attached:
        return False, "Ether not attached to glycerol"

    # Check for at least one alkyl chain (long carbon chain)
    # Look for at least 4 consecutive carbons (C-C-C-C)
    alkyl_chain = Chem.MolFromSmarts("C-C-C-C")
    if not mol.HasSubstructMatch(alkyl_chain):
        return False, "No alkyl chain detected"

    return True, "Glycerol with ether-linked alkyl chain"
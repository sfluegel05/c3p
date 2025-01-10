"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: CHEBI:XXXXX epoxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    An epoxy fatty acid is a fatty acid containing an epoxide ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (-C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Look for epoxide ring (C1OC1) with a more specific pattern
    epoxide_pattern = Chem.MolFromSmarts("[C;r3][O;r3][C;r3]")
    epoxide_matches = mol.GetSubstructMatches(epoxide_pattern)
    if not epoxide_matches:
        return False, "No epoxide ring found"

    # Ensure the epoxide ring is part of the main carbon chain
    # by checking that at least one carbon in the epoxide is connected to the rest of the molecule
    epoxide_carbons = set()
    for match in epoxide_matches:
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:
                epoxide_carbons.add(atom_idx)
    
    connected_to_rest = False
    for carbon_idx in epoxide_carbons:
        atom = mol.GetAtomWithIdx(carbon_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in epoxide_carbons:
                connected_to_rest = True
                break
        if connected_to_rest:
            break
    if not connected_to_rest:
        return False, "Epoxide ring not connected to the main carbon chain"

    # Check for a long carbon chain (at least 10 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 10:
        return False, "Carbon chain too short to be a fatty acid"

    # Check molecular weight (epoxy fatty acids typically >300 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for epoxy fatty acid"

    # Check for a typical fatty acid structure (long aliphatic chain)
    # by counting the number of rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Not enough rotatable bonds for a typical fatty acid"

    return True, "Contains a fatty acid chain with an epoxide ring"
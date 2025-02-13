"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
"""
Classifies: CHEBI:32879 arenecarbaldehyde
Any aldehyde in which the carbonyl group is attached to an aromatic moiety.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an arenecarbaldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Detect aldehydes
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if not aldehyde_matches:
        return False, "No aldehyde group found"

    # Check if there is exactly one aldehyde group
    if len(aldehyde_matches) != 1:
        return False, f"Found {len(aldehyde_matches)} aldehyde groups, expected exactly 1"

    # Get the carbonyl carbon index
    carbonyl_carbon_idx = aldehyde_matches[0][0]
    carbonyl_carbon = mol.GetAtomWithIdx(carbonyl_carbon_idx)

    # Check if the carbonyl carbon is directly attached to an aromatic ring
    if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in carbonyl_carbon.GetNeighbors()):
        return True, "Aldehyde group directly attached to an aromatic ring"

    # Check if the carbonyl carbon is connected to an aromatic moiety through a carbon-carbon bond
    aromatic_rings = [ring for ring in mol.GetRingInfo().AtomRings() if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)]
    for ring in aromatic_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6 and any(neighbor_idx == carbonyl_carbon_idx for neighbor_idx in atom.GetNeighbors()):
                return True, "Aldehyde group connected to an aromatic moiety through a carbon-carbon bond"

    # Additional checks
    if mol.GetAtomWithIdx(carbonyl_carbon_idx).GetTotalNumHs(onlyExplicit=False) != 1:
        return False, "Carbonyl group is not an aldehyde"

    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
        return False, "Molecular weight too low for an arenecarbaldehyde"

    if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 and mol.GetAtomWithIdx(idx).GetTotalNumHs(onlyExplicit=False) == 1 for idx in carbonyl_carbon.GetNeighbors()):
        return False, "Carbonyl group is part of a carboxylic acid"

    return False, "Aldehyde group not attached to an aromatic moiety"
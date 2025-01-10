"""
Classifies: CHEBI:23824 diol
"""
from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    A diol is defined as a compound containing exactly two hydroxy groups that are part of an aliphatic structure,
    not part of aromatic rings or other complex functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a simple diol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for hydroxy group pattern (-OH) on aliphatic carbon
    hydroxy_pattern = Chem.MolFromSmarts("[#6;!$(C=C)][OX2H]")  # Aliphatic carbon with attached hydroxyl
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    
    # Count the number of hydroxy groups
    num_hydroxy_groups = len(hydroxy_matches)
    
    # A diol should have exactly two -OH groups
    if num_hydroxy_groups != 2:
        return False, f"Contains {num_hydroxy_groups} hydroxy groups, need exactly 2"

    # Additional check to confirm hydroxyls are not part of aromatic or complex groups
    for match in hydroxy_matches:
        hydroxy_atom = mol.GetAtomWithIdx(match[0])
        carbon_atom = hydroxy_atom.GetNeighbors()[0]
        
        # Ensure the carbon is part of an aliphatic chain or ring
        if not carbon_atom.GetIsAromatic() and carbon_atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
            continue
        else:
            return False, "Hydroxyl groups are part of complex or aromatic groups, not simple diol"

    return True, "Contains exactly two hydroxy groups forming a simple diol"
"""
Classifies: CHEBI:24026 fatty alcohol
"""
from rdkit import Chem

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.
    A fatty alcohol is an aliphatic alcohol with a carbon chain of 3 to >27 atoms, which may be saturated/unsaturated, and branched or unbranched.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a fatty alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for aromaticity: Fatty alcohols must be purely aliphatic
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Aromatic system present, not a simple fatty alcohol"

    # Count total carbon atoms in the molecule
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 3:
        return False, f"Too few carbon atoms ({carbon_count}), need at least 3"

    # Check for sp3 hybridized carbon bound to hydroxyl groups (O-H) or primary aliphatic alcohols   
    oh_pattern = Chem.MolFromSmarts("[CX4;!$(C=[O,N])][OX2H]")  # Allow multiple hydroxyls for diols
    if not mol.HasSubstructMatch(oh_pattern):
        return False, "No aliphatic hydroxyl group found on a carbon chain"

    # Allow multiple aliphatic oxygens, reject non-carbons in excess
    non_carbon_hetero_pattern = Chem.MolFromSmarts("[!C;!H]")
    heteroatoms_outside_oxygen = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [6, 1, 8]]
    if len(heteroatoms_outside_oxygen) > 0:
        return False, f"Unexpected heteroatoms present, functionality beyond simple hydroxyl groups"

    # Allow ketones/aldehydes if not bonded to hydroxyl bearing carbons; Reliable matches with wider interpretation
    if not any(atom.GetSymbol() == 'O' and atom.GetDegree() > 1 for atom in mol.GetAtoms()):
        carbonyl_near_oh_pattern = Chem.MolFromSmarts("[CX3]=O.[OH1]")
        if mol.HasSubstructMatch(carbonyl_near_oh_pattern):
            return False, "Contains complex carbonyl groups suggesting non-primary alcohol"

    return True, "Contains carbon chain with hydroxyl group(s), consistent with fatty alcohols"
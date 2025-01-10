"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
"""
Classifies: aliphatic alcohol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    An aliphatic alcohol is an alcohol derived from an aliphatic compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for hydroxyl groups (excluding carboxylic acids)
    alcohol_pattern = Chem.MolFromSmarts("[CX4][OH1]")
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "No aliphatic hydroxyl groups found"
    
    # Count hydroxyl groups
    hydroxyl_matches = len(mol.GetSubstructMatches(alcohol_pattern))
    
    # Count carboxylic acid groups
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OH1]")
    carboxyl_matches = len(mol.GetSubstructMatches(carboxyl_pattern))
    
    # If all OH groups are from carboxylic acids, it's not an alcohol
    if hydroxyl_matches == 0:
        return False, "No alcohol groups present"
    
    # Calculate aromaticity
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    total_atoms = mol.GetNumAtoms()
    aromatic_ratio = aromatic_atoms / total_atoms if total_atoms > 0 else 0
    
    # If the molecule is primarily aromatic (>50%), reject it
    if aromatic_ratio > 0.5:
        return False, "Molecule is primarily aromatic"
    
    # Get all carbons attached to OH groups
    carbons_with_oh = []
    for match in mol.GetSubstructMatches(alcohol_pattern):
        c_idx = match[0]  # Carbon index
        carbons_with_oh.append(mol.GetAtomWithIdx(c_idx))
    
    # Check if any carbon with OH is sp3 hybridized
    has_sp3_alcohol = False
    for carbon in carbons_with_oh:
        if carbon.GetHybridization() == Chem.HybridizationType.SP3:
            has_sp3_alcohol = True
            break
            
    if not has_sp3_alcohol:
        return False, "No sp3 carbons with hydroxyl groups found"
    
    # Success case
    if hydroxyl_matches == 1:
        return True, "Contains one aliphatic hydroxyl group"
    else:
        return True, f"Contains {hydroxyl_matches} aliphatic hydroxyl groups"
"""
Classifies: CHEBI:134179 volatile organic compound
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a volatile organic compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate molecular weight
    mol_wt = Descriptors.ExactMolWt(mol)  # More accurate calculation
    
    # New criteria for volatile organic compounds
    # Allow higher molecular weight for long chains and specific functional groups
    if mol_wt > 500:
        return False, "Molecular weight exceeds typical range for volatility"

    # Check presence of aromatic rings
    aromatic_count = sum(1 for atom in mol.GetAromaticAtoms())
    
    # Check for specific volatile functional groups: alcohols, alkenes, aldehydes
    alcohol_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3](O)")))
    alkene_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3]=[CX3]")))
    aldehyde_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3H](=O)" )))
    
    # Simple rules to consider more volatility-related structures 
    if (alcohol_count > 0 or alkene_count > 0 or aldehyde_count > 0 or aromatic_count > 0) or mol_wt <= 300:
        return True, "Likely a volatile organic compound based on molecular structure and composition"
    else:
        return False, "Does not match enhanced characteristics of a volatile organic compound"
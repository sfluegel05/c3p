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
    mol_wt = Descriptors.ExactMolWt(mol)
    
    # Check molecular weight for VOC typical range
    if mol_wt > 300:
        return False, "Molecular weight exceeds typical range for volatile organic compounds"

    # Check for specific volatile functional groups
    aromatic_count = sum(1 for atom in mol.GetAromaticAtoms())
    alcohol_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3](O)")))
    alkene_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3]=[CX3]")))
    aldehyde_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3H](=O)")))
    ketone_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3](=O)[#6]")))
    halogen_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[#6][F,Cl,Br,I]")))
    ether_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3]O[CX3]")))
    
    # Additional criteria for volatility
    unsaturation_degree = alkene_count + aromatic_count
    
    # Enhanced logic for classifying VOCs
    if (
        alcohol_count > 0 or 
        aldehyde_count > 0 or 
        ketone_count > 0 or
        ether_count > 0 or
        halogen_count > 0 or
        unsaturation_degree > 0 or 
        (mol_wt <= 200 and unsaturation_degree > 0)
    ):
        return True, "Likely a volatile organic compound based on molecular structure and composition"
    else:
        return False, "Does not match enhanced characteristics of a volatile organic compound"
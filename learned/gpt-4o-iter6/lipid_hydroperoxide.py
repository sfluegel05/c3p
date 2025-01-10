"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
from rdkit import Chem

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a lipid hydroperoxide based on its SMILES string.
    A lipid hydroperoxide is any lipid carrying one or more hydroperoxy (-OOH) substituents.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipid hydroperoxide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find hydroperoxy groups (-COOH where C is sp3 with a hydrocarbon chain)
    hydroperoxy_pattern = Chem.MolFromSmarts("[C;!R][OX2H]")
    if not mol.HasSubstructMatch(hydroperoxy_pattern):
        return False, "No hydroperoxy group found"
    
    # Check for unsaturation - a typical feature in lipid molecules
    unsaturation_pattern = Chem.MolFromSmarts("C=C")
    if mol.HasSubstructMatch(unsaturation_pattern):
        unsaturation_present = True
    else:
        unsaturation_present = False
    
    # Calculate the number of carbon atoms - Typical lipid molecules have long chains
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    if carbon_count < 12:
        return False, f"Carbon chain too short for typical lipids: found {carbon_count} carbons"
    
    # Determine if reasonable lipid structure, here we assume presence of both long chains and unsaturation.
    if carbon_count >= 12 and unsaturation_present:
        return True, "Contains hydroperoxy group(s) within a lipid structure with unsaturation"
    
    return False, "Missing characteristics of a typical lipid structure despite presence of hydroperoxy group"
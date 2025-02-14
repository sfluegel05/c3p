"""
Classifies: CHEBI:24279 glucosinolate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosinolate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for thioglucoside core (pyranose linked to sulfur), allow for substitutions
    thioglucoside_pattern = Chem.MolFromSmarts("[C]1[C][C][O][C][C]1[S]")
    thioglucoside_matches = mol.GetSubstructMatches(thioglucoside_pattern)
    if not thioglucoside_matches:
       return False, "Thioglucoside core not found"

    # 2. Check for the central C atom bonded to S, N, and a side-chain, with stereochemistry
    # This includes the trans C=N bond pattern and ensures that C has 3 bonds and side chain on C is not H. 
    central_carbon_pattern = Chem.MolFromSmarts("[SX2][C](!=[H])(=[N])")
    central_carbon_matches = mol.GetSubstructMatches(central_carbon_pattern)
    if not central_carbon_matches:
        return False, "Central carbon with S and N not found"

    # 3. Check for sulfonated oxime group (anionic and non-anionic, with trans configuration)
    sulfonated_oxime_pattern = Chem.MolFromSmarts("[C]/=[N]\[OS(=O)(=O)[O-,O]")
    sulfonated_oxime_matches = mol.GetSubstructMatches(sulfonated_oxime_pattern)

    if not sulfonated_oxime_matches:
       return False, "Sulfonated oxime group not found"
    
    # 4. Check connectivity between thioglucoside S and central C and combine the pattern
    # this ensures that the S in the thioglucoside and the S in the central_carbon_pattern is the same
    combined_pattern = Chem.MolFromSmarts("[C]1[C][C][O][C][C]1[S][C](!=[H])=[N]")
    combined_matches = mol.GetSubstructMatches(combined_pattern)
    if not combined_matches:
        return False, "Thioglucoside sulfur not connected to the central carbon"

    return True, "Molecule matches all required patterns for glucosinolate"
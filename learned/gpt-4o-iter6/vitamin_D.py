"""
Classifies: CHEBI:27300 vitamin D
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D compound based on its SMILES string.
    Vitamin D compounds are defined as fat-soluble hydroxy seco-steroids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a vitamin D compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for basic vitamin D structural features
    # Seco-steroid structure (broken ring with characteristic double bonds)
    seco_steroid_pattern = Chem.MolFromSmarts("C1CCCCC1=C/C2=C1CCCC2")
    if not mol.HasSubstructMatch(seco_steroid_pattern):
        return False, "No seco-steroid backbone structure typical of vitamin D"
    
    # Check for presence of hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")  # OH group
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 1:
        return False, "Missing hydroxyl groups indicative of vitamin D"

    # Check molecular properties indicative of fat solubility
    mol_logP = rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
    if mol_logP < 3:  # Typical vitamin D has a higher logP
        return False, "Molecule not hydrophobic enough for vitamin D"
   
    # For identified compounds: typically see multiple hydroxyls
    if len(hydroxyl_matches) < 2:
        return False, f"Only {len(hydroxyl_matches)} hydroxyl groups, common vitamin D compounds have more"

    return True, "Matches typical vitamin D seco-steroid structure with hydroxyl groups"
"""
Classifies: CHEBI:76575 monoradylglycerol
"""
from rdkit import Chem

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    A monoradylglycerol is a glycerol with one acyl, alkyl, or alk-1-enyl substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Define glycerol backbone more flexibly (allowing for substitution on each hydroxyl group)
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Check for one substituent group - acyl (C=O), alkyl (C), or alk-1-enyl (C=C)
    # Ensure patterns specifically capture correct substituent connectivity
    acyl_pattern = Chem.MolFromSmarts("C(=O)O[C]")
    alkyl_pattern = Chem.MolFromSmarts("COC[C]")
    alkenyl_pattern = Chem.MolFromSmarts("COC=C[C]")

    # Count matches for each type of substituent, accounting that there should only be one substituent in total
    acyl_matches = len(mol.GetSubstructMatches(acyl_pattern))
    alkyl_matches = len(mol.GetSubstructMatches(alkyl_pattern))
    alkenyl_matches = len(mol.GetSubstructMatches(alkenyl_pattern))
    
    # Total substituent matches should be one
    if (acyl_matches + alkyl_matches + alkenyl_matches) != 1:
        return False, f"Expected 1 substituent group, found {acyl_matches + alkyl_matches + alkenyl_matches}"

    # Ensure valid lipid chain length, total carbon.
    # Here summed to ensure that a single chain gets counted as such
    glycerol_carbon_pattern = Chem.MolFromSmarts("OCC(O)CO[$([CH])]")
    c_chain_matches = mol.GetSubstructMatches(glycerol_carbon_pattern)
    
    if len(c_chain_matches) and sum([len(match) for match in c_chain_matches]) < 13:  # typical lipid length
        return False, "Chain too short for typical lipid substituent"
    
    return True, "Contains glycerol backbone with one acyl, alkyl, or alk-1-enyl substituent"
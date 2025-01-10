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

    # Check for glycerol backbone ([*]COC(CO)O[*] without specifying what attaches to each -OH)
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Count total substituents on glycerol backbone that are not part of the core
    acyl_pattern = Chem.MolFromSmarts("C(=O)O")
    alkyl_pattern = Chem.MolFromSmarts("COC")
    alkenyl_pattern = Chem.MolFromSmarts("COC=C")

    # A single substituent of any type (more flexible, focusing on connections)
    substituents = [mol.GetSubstructMatches(acyl_pattern),
                    mol.GetSubstructMatches(alkyl_pattern),
                    mol.GetSubstructMatches(alkenyl_pattern)]
    
    total_substituents = sum(len(matches) for matches in substituents)

    if total_substituents != 1:
        return False, f"Expected 1 substituent group, found {total_substituents}" 
    
    # Ensuring only one substituent is long enough to count as a separate 'radyl' group
    substituent_length_threshold = 8  # Approximate length to count as a lipid chain
    valid_substituent = any(len(matches) >= substituent_length_threshold for matches in substituents)
    if not valid_substituent:
        return False, f"No valid lipid chain substituent found"

    return True, "Contains glycerol backbone with one acyl, alkyl, or alk-1-enyl substituent"
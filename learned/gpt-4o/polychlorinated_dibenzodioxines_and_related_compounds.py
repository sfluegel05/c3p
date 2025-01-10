"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
from rdkit import Chem

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Classifies a molecule as a polychlorinated dibenzodioxin or related compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule belongs to the class, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for the relevant structures
    # General pattern for dioxins/furans (2 benzene rings connected through oxygen)
    dioxin_furan_pattern = Chem.MolFromSmarts('Oc1c2ccccc2Oc3c1ccccc3')
    
    # General pattern for biphenyl (two phenyl groups cross-coupled)
    biphenyl_pattern = Chem.MolFromSmarts('c1ccccc1-c2ccccc2')
    
    # Check for chlorination/bromination
    halogenated_pattern = Chem.MolFromSmarts('[*;#6][Cl,Br]') # Matches C-Cl or C-Br bonds.
    
    # Check for required structures
    has_dioxin_furan = mol.HasSubstructMatch(dioxin_furan_pattern)
    has_biphenyl = mol.HasSubstructMatch(biphenyl_pattern)
    
    if not (has_dioxin_furan or has_biphenyl):
        return False, "Does not contain dioxin, furan, or biphenyl structure"
    
    # Check if structure is halogenated (containing at least one C-Cl or C-Br bond)
    if not mol.HasSubstructMatch(halogenated_pattern):
        return False, "Structure is not halogenated with Cl or Br"
    
    # All checks passed
    return True, "Classified as polychlorinated dibenzodioxin or related compound"

# Example usage
example_smiles = "Clc1ccc(cc1)-c1cc(Cl)c(Cl)cc1Cl"
result, reason = is_polychlorinated_dibenzodioxines_and_related_compounds(example_smiles)
print(result, reason)
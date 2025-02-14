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
    
    # Define SMARTS patterns for typical structures:
    # Dioxins: Two benzene rings connected by two oxygen atoms.
    dioxin_pattern = Chem.MolFromSmarts('O1c2ccccc2Oc3ccccc13')
    
    # Dibenzofurans: Two benzene rings connected by one bridge oxygen.
    furan_pattern = Chem.MolFromSmarts('O1c2ccccc2-c3c1ccccc3')
    
    # Halogenated biphenyl patterns should be more specific regarding chlorine/bromine
    halogenated_biphenyl_pattern = Chem.MolFromSmarts('c1cc([Cl,Br])c([Cl,Br])c(c1)-c1c([Cl,Br])cc([Cl,Br])cc1')
    
    # Check for polychlorination/polybromination
    high_halogens_pattern = Chem.MolFromSmarts('[Cl,Br].[Cl,Br].[Cl,Br]')  # Requires at least 3 halogens
    
    # Determine the presence of the key structure motifs
    has_dioxin = mol.HasSubstructMatch(dioxin_pattern)
    has_furan = mol.HasSubstructMatch(furan_pattern)
    is_halogenated_biphenyl = mol.HasSubstructMatch(halogenated_biphenyl_pattern)
    
    # Check for high levels of halogen
    has_high_halogens = mol.HasSubstructMatch(high_halogens_pattern)
    
    # Ensure at least one significant structure type and significant halogenation
    if (has_dioxin or has_furan or is_halogenated_biphenyl) and has_high_halogens:
        return True, "Classified as polychlorinated dibenzodioxin or related compound"
    
    return False, "Does not meet classification criteria for polychlorinated dibenzodioxin or related compound"

# Example usage
example_smiles = "Clc1ccc(cc1)-c1cc(Cl)c(Cl)cc1Cl"
result, reason = is_polychlorinated_dibenzodioxines_and_related_compounds(example_smiles)
print(result, reason)
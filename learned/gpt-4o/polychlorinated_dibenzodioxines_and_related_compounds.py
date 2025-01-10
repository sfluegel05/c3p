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
    
    # Define SMARTS patterns for relevant structures
    # Dioxin pattern (basic 2-benzene rings with dioxin bridge pattern)
    dioxin_pattern = Chem.MolFromSmarts('Oc1ccc2cccc3Oc1c(c34)ccc4')
    
    # Benzofuran pattern
    benzofuran_pattern = Chem.MolFromSmarts('oc1ccc2ccccc2c1')
    
    # Biphenyl pattern
    biphenyl_pattern = Chem.MolFromSmarts('c1ccc(c2ccccc2)c1')
    
    # Check for chlorination/bromination
    halogenated_pattern = Chem.MolFromSmarts('[Cl,Br]')
    
    # Check for any benzodioxin, benzofuran or biphenyl structure
    has_dioxin = mol.HasSubstructMatch(dioxin_pattern)
    has_benzofuran = mol.HasSubstructMatch(benzofuran_pattern)
    has_biphenyl = mol.HasSubstructMatch(biphenyl_pattern)
    
    # Check for any halogen substitutions
    halogenated_atoms = mol.GetAtomWithIdx(0).GetNeighbors()
    
    if not (has_dioxin or has_benzofuran or has_biphenyl):
        return False, "Does not contain dioxin, benzofuran, or biphenyl structure"
    
    # Check if structure is halogenated (with chlorine or bromine)
    if not any(atom.GetSymbol() in ["Cl", "Br"] for atom in mol.GetAtoms()):
        return False, "Structure is not halogenated"
    
    # All checks passed
    return True, "Classified as polychlorinated dibenzodioxin or related compound"

# Example usage
example_smiles = "Clc1ccc(cc1)-c1cc(Cl)c(Cl)cc1Cl"
result, reason = is_polychlorinated_dibenzodioxines_and_related_compounds(example_smiles)
print(result, reason)
"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Classifies a molecule as a polychlorinated dibenzodioxin or a related compound based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule belongs to the class, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for chlorine atoms
    chlorine_count = sum(atom.GetSymbol() == 'Cl' for atom in mol.GetAtoms())
    if chlorine_count < 2:
        return False, "Not enough chlorine atoms"

    # SMARTS pattern for a polychlorinated biphenyl (biphenyl structure with chlorine)
    pcb_pattern = Chem.MolFromSmarts('c1ccccc1-c2ccccc2')
    
    # SMARTS pattern for dioxin/furan-like structures with chlorine substituents
    dioxin_pattern = Chem.MolFromSmarts('c1cc2Oc3ccccc3Oc2cc1')
    furan_pattern = Chem.MolFromSmarts('c1cc2oc3ccccc3c2cc1')

    # Check for presence of these patterns in the molecule
    if mol.HasSubstructMatch(pcb_pattern) or mol.HasSubstructMatch(dioxin_pattern) or mol.HasSubstructMatch(furan_pattern):
        return True, "Matches polychlorinated dibenzodioxin or related compound pattern"

    return False, "Does not match polychlorinated dibenzodioxin or related compound pattern"

# Example to test the function
smiles_example = "Clc1cc2oc3ccccc3c2cc1"
result, reason = is_polychlorinated_dibenzodioxines_and_related_compounds(smiles_example)
print(f"Result: {result}, Reason: {reason}")
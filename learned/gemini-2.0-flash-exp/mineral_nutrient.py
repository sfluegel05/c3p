"""
Classifies: CHEBI:167164 mineral nutrient
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_mineral_nutrient(smiles: str):
    """
    Determines if a molecule is a mineral nutrient based on its SMILES string.
    Mineral nutrients are inorganic compounds that provide essential elements for the body.
    They generally contain metals or other charged inorganic compounds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral nutrient, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of a metal atom
    metal_atoms = [atom.GetAtomicNum() for atom in mol.GetAtoms() if atom.GetAtomicNum() in [3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 55, 56, 57, 58, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 87, 88, 89, 90, 91, 92]]
    if not metal_atoms:
        #If no metal, check for common counterions (halides, etc) and common minerals.
        halide_pattern = Chem.MolFromSmarts("[F-,Cl-,Br-,I-]")
        sulfate_pattern = Chem.MolFromSmarts("[O-]S([O-])(=O)=O")
        phosphate_pattern = Chem.MolFromSmarts("[O-]P([O-])([O-])=O")
        carbonate_pattern = Chem.MolFromSmarts("[O-]C([O-])=O")
        nitrate_pattern = Chem.MolFromSmarts("[O-][N+]([O-])=O")
        silicate_pattern = Chem.MolFromSmarts("[O-][Si]([O-])([O-])[O-]")
        hydroxide_pattern = Chem.MolFromSmarts("[OH-]")

        if not (mol.HasSubstructMatch(halide_pattern) or mol.HasSubstructMatch(sulfate_pattern) or
                mol.HasSubstructMatch(phosphate_pattern) or mol.HasSubstructMatch(carbonate_pattern) or
                mol.HasSubstructMatch(nitrate_pattern) or mol.HasSubstructMatch(silicate_pattern) or
                mol.HasSubstructMatch(hydroxide_pattern) ):

           return False, "No metal, halide, or other common mineral found"

    # Check for carbon-hydrogen bonds, but allow for acetate, formate, or other small organic salts
    carbon_hydrogen_pattern = Chem.MolFromSmarts("[CH]")

    if mol.HasSubstructMatch(carbon_hydrogen_pattern):
        acetate_pattern = Chem.MolFromSmarts("CC([O-])=O")
        formate_pattern = Chem.MolFromSmarts("[H]C([O-])=O")

        if not (mol.HasSubstructMatch(acetate_pattern) or mol.HasSubstructMatch(formate_pattern)):
          return False, "Contains C-H bonds that are not from small salts"
    
    # Basic validation of size: Mineral Nutrients are often salts, but not large organics.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 700:
        return False, "Molecular weight too high for a mineral nutrient"

    return True, "Contains metal ions and/or typical mineral ions"
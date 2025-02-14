"""
Classifies: CHEBI:72544 flavonoids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    A flavonoid contains the C6-C3-C6 core structure, typically as a benzopyran derivative,
    and includes various subgroups such as flavones, flavanones, isoflavones, chalcones, etc.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Flavonoid substructure SMARTS patterns
    flavone_pattern = Chem.MolFromSmarts("c1cc2c(cc1)oc(=O)c(c2)c3ccccc3") # Flavone core
    flavonol_pattern = Chem.MolFromSmarts("c1cc2c(cc1)oc(=O)c(c2)c3cc(O)cc(O)c3") # Flavonol core
    flavanone_pattern = Chem.MolFromSmarts("c1cc2c(cc1)c(=O)c[cH]2c3ccccc3")  # Flavanone core
    isoflavone_pattern = Chem.MolFromSmarts("c1ccc2c(c1)c(=O)cc(c2)c3ccccc3")  # Isoflavone core
    chalcone_pattern1 = Chem.MolFromSmarts("c1ccccc1C(=O)C=Cc2ccccc2") # Chalcone
    chalcone_pattern2 = Chem.MolFromSmarts("c1ccccc1C(=O)CC=Cc2ccccc2") # Chalcone
    chalcone_pattern3 = Chem.MolFromSmarts("c1ccccc1C(=O)CCc2ccccc2") # Chalcone
    dihydrochalcone_pattern = Chem.MolFromSmarts("c1ccccc1C(=O)CCc2ccccc2") # Dihydrochalcone

    # Check if the molecule matches any of the patterns
    if (mol.HasSubstructMatch(flavone_pattern) or
        mol.HasSubstructMatch(flavonol_pattern) or
        mol.HasSubstructMatch(flavanone_pattern) or
        mol.HasSubstructMatch(isoflavone_pattern) or
        mol.HasSubstructMatch(chalcone_pattern1) or
        mol.HasSubstructMatch(chalcone_pattern2) or
        mol.HasSubstructMatch(chalcone_pattern3) or
        mol.HasSubstructMatch(dihydrochalcone_pattern)):
       return True, "Matches a flavonoid substructure pattern"
    
    return False, "Does not match any known flavonoid substructure pattern"
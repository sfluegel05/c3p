"""
Classifies: CHEBI:72544 flavonoids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    A flavonoid contains the C6-C3-C6 core structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # General Flavonoid core C6-C3-C6 pattern, allowing for substitutions
    # This pattern looks for two aromatic rings (c) connected by a 3-atom chain (~).
    # The atoms in the 3-atom chain can be C or O.
    core_pattern = Chem.MolFromSmarts("c1ccccc1~[C,O]~[C,O]~[C,O]~c2ccccc2")

    # Additional patterns for common variations - modified from the previous attempt
    # These include some specific cases that were missed or wrongly classified
    flavone_pattern = Chem.MolFromSmarts("c1cc2c(cc1)oc(=O)c(c2)c3ccccc3") # Flavone core
    flavonol_pattern = Chem.MolFromSmarts("c1cc2c(cc1)oc(=O)c(c2)c3cc(O)cc(O)c3") # Flavonol core
    flavanone_pattern = Chem.MolFromSmarts("c1cc2c(cc1)c(=O)c[cH]2c3ccccc3")  # Flavanone core
    isoflavone_pattern = Chem.MolFromSmarts("c1ccc2c(c1)c(=O)cc(c2)c3ccccc3")  # Isoflavone core
    chalcone_pattern1 = Chem.MolFromSmarts("c1ccccc1C(=O)C=Cc2ccccc2") # Chalcone
    chalcone_pattern2 = Chem.MolFromSmarts("c1ccccc1C(=O)CC=Cc2ccccc2") # Chalcone
    chalcone_pattern3 = Chem.MolFromSmarts("c1ccccc1C(=O)CCc2ccccc2") # Chalcone
    dihydrochalcone_pattern = Chem.MolFromSmarts("c1ccccc1C(=O)CCc2ccccc2") # Dihydrochalcone
    aurone_pattern = Chem.MolFromSmarts("c1cc2c(cc1)oc(c2)=Cc3ccccc3") #Aurone
    
    # Check if the molecule matches the general core OR the additional patterns
    if (mol.HasSubstructMatch(core_pattern) or
        mol.HasSubstructMatch(flavone_pattern) or
        mol.HasSubstructMatch(flavonol_pattern) or
        mol.HasSubstructMatch(flavanone_pattern) or
        mol.HasSubstructMatch(isoflavone_pattern) or
        mol.HasSubstructMatch(chalcone_pattern1) or
        mol.HasSubstructMatch(chalcone_pattern2) or
        mol.HasSubstructMatch(chalcone_pattern3) or
        mol.HasSubstructMatch(dihydrochalcone_pattern) or
        mol.HasSubstructMatch(aurone_pattern)
        ):
        return True, "Matches a flavonoid substructure pattern"
    
    return False, "Does not match any known flavonoid substructure pattern"
"""
Classifies: CHEBI:192499 anthoxanthin
"""
"""
Classifies: CHEBI:26580 anthoxanthin
Anthoxanthins are a type of flavonoid pigments in plants. They are water-soluble pigments which range in color from white or colorless to a creamy to yellow, often on petals of flowers.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthoxanthin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for flavonoid backbone
    flavonoid_patterns = [
        Chem.MolFromSmarts("c1c(oc2ccccc2)cc2ccccc12"),  # Flavone
        Chem.MolFromSmarts("c1c(oc2ccccc2)cc2cccc(O)c12"),  # Flavonol
        Chem.MolFromSmarts("c1c(oc2ccccc2)ccc1O"),  # Flavanone
        Chem.MolFromSmarts("c1c(oc2ccccc2)cc(O)c1O"),  # Flavan-3-ol
        Chem.MolFromSmarts("c1c(oc2ccccc2)ccc1"),  # Flavene
        Chem.MolFromSmarts("c1c(oc2ccccc2)ccc1O"),  # Flavanone
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in flavonoid_patterns):
        return False, "No flavonoid backbone found"

    # Check for anthoxanthin-specific substituents
    anthoxanthin_patterns = [
        Chem.MolFromSmarts("O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"),  # Glucose
        Chem.MolFromSmarts("O[C@@H](CO)[C@H](O)[C@H](O)[C@H]1O"),  # Rhamnoside
        Chem.MolFromSmarts("OS(O)(=O)=O"),  # Sulfate
        Chem.MolFromSmarts("OC"),  # Methoxy
        Chem.MolFromSmarts("O"),  # Hydroxy
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in anthoxanthin_patterns):
        return False, "No anthoxanthin-specific substituents found"

    # Check molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 800:
        return False, "Molecular weight outside typical range for anthoxanthins"

    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable > 10:
        return False, "Too many rotatable bonds for an anthoxanthin"

    logp = rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
    if logp < -2 or logp > 5:
        return False, "logP outside typical range for anthoxanthins"

    return True, "Meets structural and molecular property criteria for anthoxanthins"
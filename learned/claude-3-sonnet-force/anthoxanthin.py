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
    flavonoid_pattern = Chem.MolFromSmarts("c1c(oc2ccccc2)cc2ccccc12")  # Flavone
    if not mol.HasSubstructMatch(flavonoid_pattern):
        return False, "No flavone backbone found"

    # Check for anthoxanthin-specific substituents and patterns
    anthoxanthin_patterns = [
        Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@H]1O"),  # Glucose
        Chem.MolFromSmarts("O[C@@H](C)[C@H](O)[C@@H](O)[C@H]1O"),  # Rhamnoside
        Chem.MolFromSmarts("OS(O)(=O)=O"),  # Sulfate
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in anthoxanthin_patterns):
        return False, "No anthoxanthin-specific substituents found"

    # Check for specific substitution patterns
    anthoxanthin_subst_patterns = [
        Chem.MolFromSmarts("c1cc(O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)cc(O)c1"),  # 3-O-glycosylated
        Chem.MolFromSmarts("c1cc(OS(O)(=O)=O)cc(O)c1"),  # 3-O-sulfated
        Chem.MolFromSmarts("c1cc(O)c(OS(O)(=O)=O)c(O)c1"),  # 7-O-sulfated
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in anthoxanthin_subst_patterns):
        return False, "Substitution pattern not characteristic of anthoxanthins"

    # Check molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 800:
        return False, "Molecular weight outside typical range for anthoxanthins"

    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable > 10:
        return False, "Too many rotatable bonds for an anthoxanthin"

    logp = rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
    if logp < -2 or logp > 5:
        return False, "logP outside typical range for anthoxanthins"

    return True, "Meets structural and molecular property criteria for anthoxanthins"
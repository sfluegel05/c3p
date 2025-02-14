"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
"""
Classifies: CHEBI:62543 beta-D-glucosiduronate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_beta_D_glucosiduronate(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronate based on its SMILES string.
    A beta-D-glucosiduronate is a carbohydrate acid derivative anion obtained by
    deprotonation of the carboxy group of any beta-D-glucosiduronic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucosiduronate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glucuronide substructure pattern ([O-]C(=O)OC1OC(C(O)C(O)C1O)*)
    glucuronide_pattern = Chem.MolFromSmarts("[O-]C(=O)OC1OC(C(O)C(O)C1O)*")
    if not mol.HasSubstructMatch(glucuronide_pattern):
        return False, "No glucuronide substructure found"

    # Check for carboxylate group
    carboxylate_pattern = Chem.MolFromSmarts("[O-]C(=O)")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group found"

    # Check for sugar-like moiety
    sugar_pattern = Chem.MolFromSmarts("OC1OC(C(O)C(O)C1O)*")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar-like moiety found"

    # Additional checks
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 1000:
        return False, "Molecular weight out of typical range for beta-D-glucosiduronates"

    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Too few rotatable bonds for a typical beta-D-glucosiduronate"

    return True, "Contains a glucuronide substructure with a carboxylate group and a sugar-like moiety"
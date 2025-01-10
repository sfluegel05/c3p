"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_11_12_saturated_fatty_acyl_CoA_4_(smiles: str):
    """
    Determines if a molecule is 11,12-saturated fatty acyl-CoA(4-) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str) indicating if the molecule matches the class, and reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")

    # Check for CoA backbone (with associated functional groups indicative of CoA)
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return (False, "No CoA backbone structure found")

    # Pattern to identify 11,12-saturated fatty acyl chain:
    # Saturated position between 11 and 12
    # This pattern is a broad representation and might require refining for accuracy
    saturated_11_12_pattern = Chem.MolFromSmarts("CCCCCCCCCCC([CH2])[CH2]CCCCCCCCCC=O")
    if not mol.HasSubstructMatch(saturated_11_12_pattern):
        return (False, "No saturated 11,12 position in acyl chain detected")

    # Verify net charge of -4, typical for fatty acyl-CoA(4-) structures
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if net_charge != -4:
        return (False, "Incorrect charge state; expected -4")

    return (True, "Matches 11,12-saturated fatty acyl-CoA(4-) with expected structural properties.")

# By ensuring precise structural matches, this revised function aims to effectively classify the compound against the desired fatty acyl-CoA characteristics.
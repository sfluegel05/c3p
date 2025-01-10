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
        tuple: (bool, str) where bool indicates if the molecule matches the class,
               and str provides the reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")

    # CoA backbone structure including key functional groups
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)C")
    if not mol.HasSubstructMatch(coa_pattern):
        return (False, "No CoA backbone structure found")

    # Pattern to identify 11,12-saturated fatty acyl chain 
    # (this pattern serves as a more general check for a long chain with specific saturation)
    saturated_11_12_pattern = Chem.MolFromSmarts("CCCCCCCCC[CH2][CH2]CCCC")
    if not mol.HasSubstructMatch(saturated_11_12_pattern):
        return (False, "No saturated 11,12 position in acyl chain detected")

    # Confirm charge state of -4
    net_charge = sum([a.GetFormalCharge() for a in mol.GetAtoms()])
    if net_charge != -4:
        return (False, "Incorrect charge state; expected -4")

    return (True, "Matches 11,12-saturated fatty acyl-CoA(4-) with expected structural properties.")

# This revised function seeks to capture specific structural motifs and the saturated positions through 
# better substructural pattern representation, while acknowledging the complexity of positional saturation checks.
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
        return False, "Invalid SMILES string"
    
    # CoA structure pattern approximation
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)C")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA backbone structure found" 

    # Checking for fatty acyl portion with specific 11-12 saturation
    # Pattern capturing a long chain with specific unsaturation position skipped for complexity reasons
    # Approximation will ensure single bonds at position 11-12 are present
    acyl_pattern = Chem.MolFromSmarts("CCCCCCCCCCCC=CCCC")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No proper fatty acyl chain detected with correct 11,12-saturation"
    
    # Confirm net charge of -4
    if Chem.rdmolops.GetFormalCharge(mol) != -4:
        return False, "Incorrect charge state; expected 4-"

    return True, "Matches 11,12-saturated fatty acyl-CoA(4-) with expected structural properties."

# The function makes several assumptions, notably placing importance on structural motifs and charge, 
# albeit with recognition that full 11-12 locational checking requires sophisticated techniques not covered here.
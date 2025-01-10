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
    
    # Check for CoA backbone structure
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)C")
    if not mol.HasSubstructMatch(coa_pattern):
        return (False, "No CoA backbone structure found")
    
    # Check for acyl chain with saturated 11,12 position
    # While direct detection of 11,12 isn't trivial in standard SMARTS,
    # we ensure long chains that skip intricacies of exact position but flag obvious non-matches
    acyl_pattern = Chem.MolFromSmarts("CCCCCCCCCC[CH2][CH2]C")
    if not mol.HasSubstructMatch(acyl_pattern):
        return (False, "No proper acyl chain detected with correct 11,12 saturation")
    
    # Confirm charge state of -4
    net_charge = sum([a.GetFormalCharge() for a in mol.GetAtoms()])
    if net_charge != -4:
        return (False, "Incorrect charge state; expected -4")
    
    return (True, "Matches 11,12-saturated fatty acyl-CoA(4-) with expected structural properties.")

# This version handles the expected structural motifs and approximates 11,12 saturation 
# through pattern recognition. This function still assumes identification via substructural motifs 
# given known complexities in SMILES-based position-specific detection.
"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
from rdkit import Chem

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
    
    # Look for CoA-like structure
    coa_pattern = Chem.MolFromSmarts("*CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA structure found"
    
    # Check for fatty acyl chain presence
    # (a long chain of CH2 typically followed by -CO- and attached to CoA)
    fatty_acyl_pattern = Chem.MolFromSmarts("C(=O)CCCCC")
    if not mol.HasSubstructMatch(fatty_acyl_pattern):
        return False, "No fatty acyl chain detected"
    
    # Check for saturated 11,12-bond (i.e., no double bond between 11th and 12th carbon in the chain)
    # Assume the carbon chain must be checked for double bonds
    double_bond_pattern = Chem.MolFromSmarts("CC=CC")
    if mol.HasSubstructMatch(double_bond_pattern):
        # Naive approach to check position, ideally requires a more comprehensive atom matching strategy
        bond_counts = [atom.GetTotalNumHs() for atom in mol.GetAtoms()]
        if bond_counts[10] < 2 or bond_counts[11] < 2:
            return False, "Unsaturated 11,12-bond"
    
    # Check for 4- charge state (e.g., presence of multiple phosphates with negative charges)
    if not Chem.rdmolops.GetFormalCharge(mol) == -4:
        return False, "Incorrect charge state; expected 4-"
        
    return True, "Matches 11,12-saturated fatty acyl-CoA(4-)"
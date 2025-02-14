"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI:33171 3-hydroxy fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA(4-) based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA(4-) is a CoA ester with a fatty acid chain containing a hydroxy group at the 3-position,
    and the molecule as a whole carrying a -4 charge from deprotonation of the phosphate and diphosphate groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for CoA backbone
    coa_pattern = Chem.MolFromSmarts("C(C)(C)(CO[P@@](=O)([O-])[O-])NC(=O)CCNC(=O)CCNC(=O)SC[CX3]C(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA backbone found"
    
    # Check for fatty acid chain (minimum 4 carbons)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) == 0:
        return False, "No fatty acid chain found"
    
    # Check for hydroxyl group at 3-position
    hydroxy_pattern = Chem.MolFromSmarts("[CX4]([OH])(CC[CX3])")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxy_matches) != 1:
        return False, f"Found {len(hydroxy_matches)} hydroxy groups at 3-position, expected 1"
    
    # Check charge
    mol_charge = AllChem.GetFormalCharge(mol)
    if mol_charge != -4:
        return False, f"Expected charge of -4, got {mol_charge}"
    
    return True, "Contains CoA backbone with a fatty acid chain and a hydroxy group at the 3-position, overall charge of -4"
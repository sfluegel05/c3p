"""
Classifies: CHEBI:57347 3-oxo-fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA(4-) based on its SMILES string.
    A 3-oxo-fatty acyl-CoA(4-) is an acyl-CoA(4-) arising from deprotonation of the phosphate and diphosphate groups of any 3-oxo-fatty acyl-CoA.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a CoA moiety
    # CoA moiety includes a pantetheine group, a phosphate group, and an adenosine group
    # Using a more general pattern to capture variations
    coa_pattern = Chem.MolFromSmarts("[O-]P([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"

    # Check for the 3-oxo-fatty acid chain (a ketone at the 3rd position from the CoA attachment)
    # The ketone should be at the 3rd position from the thioester bond (S-C=O)
    oxo_fatty_acid_pattern = Chem.MolFromSmarts("[CX3](=O)CC(=O)S")
    if not mol.HasSubstructMatch(oxo_fatty_acid_pattern):
        return False, "No 3-oxo-fatty acid chain found"

    # Check for the 4- charge (deprotonated phosphate and diphosphate groups)
    # The CoA moiety should have 4 negative charges, but we relax this to at least 4 negative charges
    charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if charge > -4:
        return False, f"Expected at least 4 negative charges, found {charge}"

    # Check for the presence of a long carbon chain (fatty acid)
    # The chain should have at least 8 carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 8:
        return False, "Carbon chain too short to be a fatty acid"

    return True, "Contains 3-oxo-fatty acid chain, CoA moiety, and at least 4 negative charges"
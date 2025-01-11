"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
"""
Classifies: CHEBI:76224 3-oxo-fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.
    A 3-oxo-fatty acyl-CoA is an oxo fatty acyl-CoA that results from the formal condensation
    of the thiol group of coenzyme A with the carboxy group of any 3-oxo-fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the CoA moiety pattern
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"

    # Define the 3-oxo-fatty acid pattern (ketone at the 3rd position)
    oxo_fatty_acid_pattern = Chem.MolFromSmarts("CC(=O)CC(=O)S")
    oxo_fatty_acid_matches = mol.GetSubstructMatches(oxo_fatty_acid_pattern)
    if len(oxo_fatty_acid_matches) == 0:
        return False, "No 3-oxo-fatty acid chain found"

    # Check if the 3-oxo-fatty acid is thioesterified with the CoA thiol group
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2][CX4]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) == 0:
        return False, "No thioester bond between CoA and 3-oxo-fatty acid"

    # Count the number of carbons in the fatty acid chain
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 6:
        return False, "Fatty acid chain too short"

    # Check molecular weight - 3-oxo-fatty acyl-CoA typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for 3-oxo-fatty acyl-CoA"

    return True, "Contains CoA moiety thioesterified with a 3-oxo-fatty acid chain"
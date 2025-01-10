"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
"""
Classifies: CHEBI:28494 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA is a CoA ester of a 3-hydroxy fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for CoA moiety pattern (thiol group attached to adenine, ribose, and phosphate groups)
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"

    # Look for 3-hydroxy fatty acid chain (long carbon chain with hydroxyl at third position)
    hydroxy_fatty_acid_pattern = Chem.MolFromSmarts("[CX4][CX4][CX4]([OH])[CX4][CX4][CX4](=O)S")
    if not mol.HasSubstructMatch(hydroxy_fatty_acid_pattern):
        return False, "No 3-hydroxy fatty acid chain found"

    # Check for ester bond between CoA and fatty acid
    ester_bond_pattern = Chem.MolFromSmarts("[CX4][CX4][CX4]([OH])[CX4][CX4][CX4](=O)S")
    if not mol.HasSubstructMatch(ester_bond_pattern):
        return False, "No ester bond between CoA and fatty acid"

    # Count carbons in the fatty acid chain
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, "Fatty acid chain too short"

    # Count oxygens to ensure the presence of hydroxyl and carboxyl groups
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:
        return False, "Insufficient oxygen atoms for 3-hydroxy fatty acyl-CoA"

    return True, "Contains CoA moiety with 3-hydroxy fatty acid chain attached via ester bond"
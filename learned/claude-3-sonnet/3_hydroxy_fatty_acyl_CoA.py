"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acyl_CoA(smiles):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3-hydroxy fatty acyl-CoA, False otherwise
        str: Reason for the classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA moiety
    coa_pattern = Chem.MolFromSmarts("C(C)(COP(OP(OC[C@H]1O[C@@H]([C@H](O)[C@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)(O)=O)(O)=O)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing CoA moiety"

    # Check for thioester bond
    thioester_pattern = Chem.MolFromSmarts("CCS(=O)")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Missing thioester bond"

    # Check for 3-hydroxy group
    hydroxy_pattern = Chem.MolFromSmarts("[C@@H](O)CC")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if not hydroxy_matches:
        return False, "Missing 3-hydroxy group"

    # Check for fatty acid chain
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
        return False, "Missing fatty acid chain"

    # Count rotatable bonds to check chain length
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Fatty acid chain too short"

    return True, "Contains CoA moiety, 3-hydroxy group, and a fatty acid chain"
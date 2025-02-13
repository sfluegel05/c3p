"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acyl_CoA(smiles):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA based on its SMILES string.
    A monounsaturated fatty acyl-CoA is an unsaturated fatty acyl-CoA with one carbon-carbon double bond in the fatty acyl chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a monounsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA substructure
    coa_pattern = Chem.MolFromSmarts("C(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing CoA substructure"

    # Check for monounsaturated fatty acyl chain
    chain_patterns = [
        Chem.MolFromSmarts("[CX4,CX3]~[CX3]=[CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]"),  # Linear aliphatic chain with double bond
        Chem.MolFromSmarts("[CX4,CX3]~[CX3]=[CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]"),  # Longer linear aliphatic chain with double bond
        Chem.MolFromSmarts("[CX4,CX3]~[CX3]=[CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]"),  # Even longer linear aliphatic chain with double bond
        Chem.MolFromSmarts("[CX4,CX3]~[CX3]=[CX3]~[CX4,CX3]~[CX3](=O)"),  # Chain with ketone group
        Chem.MolFromSmarts("[CX4,CX3]~[CX3]=[CX3]~[CX4,CX3]~[OX2]"),  # Chain with hydroxyl group
        # Add more patterns as needed to cover additional functional groups or arrangements
    ]
    chain_match = False
    for pattern in chain_patterns:
        if mol.HasSubstructMatch(pattern):
            chain_match = True
            break
    if not chain_match:
        return False, "Missing monounsaturated fatty acyl chain"

    # Check for appropriate number of rotatable bonds (proxy for chain length)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Fatty acyl chain too short"

    # Check for appropriate element counts
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count < 10 or o_count < 6:
        return False, "Insufficient carbon or oxygen atoms for fatty acyl-CoA"

    return True, "Molecule is a monounsaturated fatty acyl-CoA"
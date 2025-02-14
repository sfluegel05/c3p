"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
"""
Classifies: Long-chain fatty acyl-CoA
A fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A 
with the carboxy group of any long-chain (C13 to C22) fatty acid.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for CoA backbone pattern
    coa_pattern = Chem.MolFromSmarts("C1OC(COP(=O)(O)OP(=O)(O)OCC2OC(N3C=NC4=C3N=CN=C4N)C(O)C2OP(=O)(O)O)C(O)C1OP(=O)(O)O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA backbone found"

    # Look for fatty acyl chain pattern
    fatty_acyl_patterns = [
        Chem.MolFromSmarts("[CX4]([CX4])([CX4])[CX4]"),  # Branched chain
        Chem.MolFromSmarts("[CX4H3]"),  # Terminal methyl group
        Chem.MolFromSmarts("[CX4H2]"),  # Internal methylenes
        Chem.MolFromSmarts("[CX3H1]"),  # Internal methines
        Chem.MolFromSmarts("[CX3H0]"),  # Quaternary carbons
        Chem.MolFromSmarts("[CX2H2]"),  # Alkene carbons
        Chem.MolFromSmarts("[CX1H1]"),  # Alkyne carbons
        Chem.MolFromSmarts("[OX2H1]"),  # Hydroxyl groups
        Chem.MolFromSmarts("[OX1H0]"),  # Keto groups
        Chem.MolFromSmarts("[SX2]")  # Thioester linkage
    ]
    fatty_acyl_matches = [mol.GetSubstructMatches(pattern) for pattern in fatty_acyl_patterns]
    fatty_acyl_atoms = set().union(*fatty_acyl_matches)

    # Check for fatty acyl chain length
    n_carbons = sum(1 for atom in fatty_acyl_atoms if atom.GetAtomicNum() == 6)
    if not (13 <= n_carbons <= 22):
        return False, f"Fatty acyl chain too short (C{n_carbons})"

    # Check for thioester linkage
    thioester_match = mol.GetSubstructMatch(Chem.MolFromSmarts("[SX2]"))
    if not thioester_match:
        return False, "No thioester linkage found"

    # Check for rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Fatty acyl chain too rigid"

    return True, "Contains long-chain fatty acyl group linked to CoA via a thioester bond"
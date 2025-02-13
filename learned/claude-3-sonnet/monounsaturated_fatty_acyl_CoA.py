"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
"""
Classifies: CHEBI:49823 monounsaturated fatty acyl-CoA

A monounsaturated fatty acyl-CoA is any unsaturated fatty acyl-CoA in which the fatty acyl chain
contains one carbon-carbon double bond.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_monounsaturated_fatty_acyl_CoA(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a monounsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for CoA substructure
    coa_pattern = Chem.MolFromSmarts("C(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing CoA substructure"

    # Look for monounsaturated fatty acyl chain (at least 6 carbon atoms, linear, aliphatic, with one double bond)
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX3]=[CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if len(chain_matches) < 1:
        return False, "No monounsaturated fatty acyl chain found"

    # Count rotatable bonds to verify long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
        return False, "Fatty acyl chain too short"

    # Count carbon and oxygen atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if c_count < 20:
        return False, "Too few carbon atoms for a monounsaturated fatty acyl-CoA"
    if o_count != 6:
        return False, "Expected 6 oxygen atoms for a monounsaturated fatty acyl-CoA"

    # Exclude molecules with ring structures
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Monounsaturated fatty acyl-CoAs should not contain ring structures"

    return True, "Contains a CoA group attached to a monounsaturated fatty acyl chain"
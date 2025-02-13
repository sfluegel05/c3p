"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
"""
Classifies: CHEBI:138163 3-oxo-fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_fatty_acyl_CoA(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.
    A 3-oxo-fatty acyl-CoA is defined as an oxo fatty acyl-CoA that results from the formal condensation
    of the thiol group of coenzyme A with the carboxy group of any 3-oxo-fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple[bool, str]: A tuple containing:
            bool: True if the molecule is a 3-oxo-fatty acyl-CoA, False otherwise
            str: Reason for the classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA substructure
    coa_pattern = Chem.MolFromSmarts("CNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing CoA substructure"

    # Check for 3-oxo group on fatty acid chain
    oxo_pattern = Chem.MolFromSmarts("CC(=O)")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if not oxo_matches:
        return False, "No 3-oxo group found on fatty acid chain"

    # Check for specific double bond patterns and stereochemistry
    double_bond_patterns = [
        Chem.MolFromSmarts(r"/C=C/"),     # cis double bond
        Chem.MolFromSmarts(r"/C=C\"),     # trans double bond
        Chem.MolFromSmarts(r"C=C")        # non-stereoisomeric double bond
    ]
    double_bond_matches = []
    for pattern in double_bond_patterns:
        double_bond_matches.extend(mol.GetSubstructMatches(pattern))

    # Check molecular weight and atom counts
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Classify based on various criteria
    if mol_wt >= 500 and c_count >= 20 and o_count >= 6:
        if double_bond_matches:
            reason = f"Contains CoA substructure, 3-oxo group, {len(double_bond_matches)} double bonds"
        else:
            reason = "Contains CoA substructure, 3-oxo group, and a saturated fatty acid chain"
        return True, reason
    else:
        return False, "Molecule does not meet the criteria for a 3-oxo-fatty acyl-CoA"
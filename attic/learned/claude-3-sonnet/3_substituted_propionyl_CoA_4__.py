"""
Classifies: CHEBI:65111 3-substituted propionyl-CoA(4-)
"""
"""
Classifies: 3-substituted propionyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_substituted_propionyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-substituted propionyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-substituted propionyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly 4 negative charges
    charge_pattern = Chem.MolFromSmarts("[O-]")
    charge_matches = mol.GetSubstructMatches(charge_pattern)
    if len(charge_matches) != 4:
        return False, f"Found {len(charge_matches)} negative charges, need exactly 4"

    # Check for adenine base (more flexible pattern)
    adenine_pattern = Chem.MolFromSmarts("c1ncnc2[nH0]cnc12")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "No adenine base found"

    # Check for ribose (simplified pattern)
    ribose_pattern = Chem.MolFromSmarts("OC1C(O)C(CO)OC1")
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "No ribose moiety found"

    # Check for phosphate groups (more flexible pattern)
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)([O-])O")
    phosphate_matches = len(mol.GetSubstructMatches(phosphate_pattern))
    if phosphate_matches < 3:
        return False, f"Insufficient phosphate groups found: {phosphate_matches}"

    # Check for pantetheine arm (more flexible pattern)
    pantetheine_pattern = Chem.MolFromSmarts("NCCC(=O)NCCSC(=O)")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "No pantetheine arm found"

    # Check for thioester linkage
    thioester_pattern = Chem.MolFromSmarts("SC(=O)C")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Verify presence of acyl chain with at least 3 carbons
    acyl_chain = Chem.MolFromSmarts("SC(=O)CC")
    if not mol.HasSubstructMatch(acyl_chain):
        return False, "Acyl chain too short or missing"

    # Count key atoms to verify overall composition
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)

    # Basic sanity checks on atom counts
    if p_count != 3:
        return False, f"Must have exactly 3 phosphorus atoms, found {p_count}"
    if s_count != 1:
        return False, f"Must have exactly 1 sulfur atom, found {s_count}"
    if n_count < 5:
        return False, f"Insufficient nitrogen atoms for CoA structure, found {n_count}"
    if c_count < 20:
        return False, f"Insufficient carbon atoms for CoA structure, found {c_count}"

    # All criteria met
    return True, "Contains CoA(4-) moiety with appropriate thioester linkage and substituted acyl chain"
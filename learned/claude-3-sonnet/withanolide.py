"""
Classifies: CHEBI:74716 withanolide
"""
"""
Classifies: withanolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    A withanolide is a C28 steroid lactone with a modified side chain forming 
    a lactone ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a withanolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons - withanolides should have 28 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (26 <= c_count <= 30):  # Allow some flexibility for derivatives
        return False, f"Carbon count {c_count} not consistent with withanolide structure (expected ~28)"

    # Check for steroid core (four fused rings)
    steroid_core = Chem.MolFromSmarts("[C]1[C][C]2[C][C][C]3[C][C][C]4[C][C][C][C]4[C]3[C]2[C]1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Check for lactone ring (usually δ-lactone)
    lactone_pattern = Chem.MolFromSmarts("O=C1OCC[C]1")  # Basic lactone pattern
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring found"

    # Check for typical oxygen count (usually 4-7 oxygens including lactone)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if not (3 <= o_count <= 12):  # Allow range for various derivatives
        return False, f"Oxygen count {o_count} not typical for withanolides"

    # Check for reasonable molecular weight range
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (400 <= mol_wt <= 800):  # Typical range for withanolides
        return False, f"Molecular weight {mol_wt} outside typical range for withanolides"

    # Check for ring count (should have at least 5 rings - 4 from steroid core + lactone)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 5:
        return False, f"Ring count {ring_count} too low for withanolide structure"

    # Check for specific withanolide substructure pattern
    # This pattern looks for the characteristic steroid core with lactone side chain
    withanolide_pattern = Chem.MolFromSmarts("[#6]1~[#6]~[#6]2~[#6]~[#6]~[#6]3~[#6]~[#6]~[#6]4~[#6]~[#6]~[#6]~[#6]4~[#6]3~[#6]2~[#6]1")
    if not mol.HasSubstructMatch(withanolide_pattern):
        return False, "Does not match characteristic withanolide structure"

    # If all checks pass, it's likely a withanolide
    return True, "Matches withanolide structure with steroid core and lactone ring"
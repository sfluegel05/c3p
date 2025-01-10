"""
Classifies: CHEBI:26125 phytosterols
"""
"""
Classifies: CHEBI:26125 phytosterols
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    Phytosterols are sterols similar to cholesterol which occur in plants and vary 
    only in carbon side chains and/or presence or absence of double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phytosterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More flexible steroid core pattern that accounts for variations
    # This pattern matches the basic 6-6-6-5 ring system with flexible bond types
    steroid_core = Chem.MolFromSmarts("C1C[C,=C]C2[C,=C]C[C,=C]C3[C,=C][C,=C][C,=C]4C[C,=C][C,=C][C@]4(C)[C@]3(C)[C@]2(C)C1")
    
    # Alternative core pattern for 9,19-cyclo variants
    cyclo_core = Chem.MolFromSmarts("C12C[C,=C]C3[C,=C]C[C,=C]C4[C,=C][C,=C][C,=C](C4)[C@]3(C)[C@]1(C)CC2")
    
    if not (mol.HasSubstructMatch(steroid_core) or mol.HasSubstructMatch(cyclo_core)):
        return False, "No steroid core structure found"

    # Check for 3β-hydroxyl group specifically
    beta_hydroxyl = Chem.MolFromSmarts("[H][C@@]1([C,=C][C,=C][C,=C]2)C[C@@H](O)CC[C@]12C")
    if not mol.HasSubstructMatch(beta_hydroxyl):
        return False, "No 3β-hydroxyl group found"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Phytosterols typically have 27-30 carbons
    if c_count < 27 or c_count > 35:
        return False, f"Carbon count ({c_count}) outside typical range for phytosterols (27-35)"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350 or mol_wt > 650:  # Increased upper limit to account for glycosides
        return False, "Molecular weight outside typical range for phytosterols"

    # Check for characteristic side chain at C17
    side_chain = Chem.MolFromSmarts("[CH2,CH3][CH,C]([CH3,CH2])[CH2,CH][CH2,CH][C,CH]")
    if not mol.HasSubstructMatch(side_chain):
        return False, "No characteristic phytosterol side chain found"

    # Check for characteristics that would exclude it from being a phytosterol
    allowed_elements = {1, 6, 8} # H, C, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_elements:
            return False, "Contains elements other than C, H, and O"

    # Count rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings"

    return True, "Contains steroid core with 3β-hydroxyl group and characteristic phytosterol side chain"
"""
Classifies: CHEBI:134179 volatile organic compound
"""
"""
Classifies: CHEBI:33262 volatile organic compound (VOC)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound based on its SMILES string.
    VOCs are organic compounds with boiling point ≤ 250°C at standard pressure (101.3 kPa).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a VOC, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if organic (contains carbon)
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Not an organic compound (no carbon atoms)"

    # Calculate molecular descriptors
    mol_wt = Descriptors.ExactMolWt(mol)
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    tpsa = Descriptors.TPSA(mol)
    log_p = Descriptors.MolLogP(mol)
    
    # Count atoms and rings
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    num_heteroatoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6])
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    
    # Structural patterns
    alcohol_pattern = Chem.MolFromSmarts("[CH2,CH3,CH1][OH1]")
    alkene_pattern = Chem.MolFromSmarts("[C]=[C]")
    alkyne_pattern = Chem.MolFromSmarts("[C]#[C]")
    small_ring_pattern = Chem.MolFromSmarts("[r3,r4,r5,r6]")
    complex_ring_pattern = Chem.MolFromSmarts("[r7,r8,r9,r10,r11,r12]")
    
    has_alcohol = mol.HasSubstructMatch(alcohol_pattern)
    has_alkene = mol.HasSubstructMatch(alkene_pattern)
    has_alkyne = mol.HasSubstructMatch(alkyne_pattern)
    has_small_ring = mol.HasSubstructMatch(small_ring_pattern)
    has_complex_ring = mol.HasSubstructMatch(complex_ring_pattern)

    # Automatic rejections
    if mol_wt > 350 and not has_alcohol:  # Higher limit for alcohols
        return False, f"Molecular weight ({mol_wt:.1f}) too high for VOC"
    
    if tpsa > 90:  # Reduced from 100
        return False, f"Too polar (TPSA={tpsa:.1f}) for VOC"
        
    if num_rings > 4 and not has_alcohol:
        return False, "Too many ring systems for VOC"
        
    if has_complex_ring and num_heteroatoms > 3:
        return False, "Complex ring system with multiple heteroatoms"
        
    if num_aromatic_rings > 3:
        return False, "Too many aromatic rings for VOC"

    # Special cases - always consider as VOCs
    if any([
        (mol_wt < 150),  # Very small molecules
        (mol_wt < 200 and num_rings == 1),  # Small monocyclic compounds
        (has_alcohol and num_carbons <= 20),  # Medium-chain alcohols
        (has_alkene and mol_wt < 250 and num_rings <= 2),  # Small alkenes
        (has_alkyne and mol_wt < 200),  # Small alkynes
    ]):
        reasons = []
        if mol_wt < 150:
            reasons.append(f"very low molecular weight ({mol_wt:.1f})")
        if mol_wt < 200 and num_rings == 1:
            reasons.append("small monocyclic compound")
        if has_alcohol and num_carbons <= 20:
            reasons.append("medium-chain alcohol")
        if has_alkene and mol_wt < 250:
            reasons.append("small alkene")
        if has_alkyne and mol_wt < 200:
            reasons.append("small alkyne")
        return True, "VOC due to: " + ", ".join(reasons)

    # Long-chain alcohol special case
    if has_alcohol and num_carbons <= 30:  # Increased from previous limit
        if rotatable_bonds <= 25:  # Increased from 20
            return True, "Long-chain alcohol within VOC parameters"
    
    # Combined property score for other cases
    volatility_score = (
        -0.004 * mol_wt +  # Reduced weight factor
        -0.2 * rotatable_bonds +  # Reduced penalty for flexibility
        -0.05 * tpsa +  # Reduced polarity penalty
        -0.1 * log_p +  # Reduced hydrophobicity penalty
        -1.0 * num_rings +  # Penalty for ring systems
        -0.5 * num_heteroatoms +  # Penalty for heteroatoms
        (2 if has_alcohol else 0) +  # Bonus for alcohol groups
        (1 if has_alkene else 0) +  # Bonus for alkenes
        (1 if has_alkyne else 0) +  # Bonus for alkynes
        (0.5 if has_small_ring else 0)  # Reduced bonus for rings
    )
    
    if volatility_score > -8:  # Adjusted threshold
        return True, "VOC based on combined physicochemical properties"
        
    return False, "Properties suggest non-volatile compound"
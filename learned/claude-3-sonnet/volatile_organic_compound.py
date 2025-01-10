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
    
    # Count specific atoms and features
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    num_heavy_atoms = mol.GetNumHeavyAtoms()
    
    # Rules for non-volatile compounds
    
    # Very large molecules are unlikely to be volatile
    if mol_wt > 400:
        return False, f"Molecular weight ({mol_wt:.1f}) too high for VOC"
    
    # Many rotatable bonds indicate flexibility and higher boiling point
    if rotatable_bonds > 20:
        return False, f"Too many rotatable bonds ({rotatable_bonds}) for VOC"
    
    # High polarity (TPSA) generally means higher boiling point
    if tpsa > 100:
        return False, f"Too polar (TPSA={tpsa:.1f}) for VOC"

    # Extremely hydrophobic compounds (high logP) tend to have high boiling points
    if log_p > 8:
        return False, f"Too hydrophobic (logP={log_p:.1f}) for VOC"

    # Common VOC structural patterns
    alcohol_pattern = Chem.MolFromSmarts("[CH2,CH3,CH1][OH1]")
    alkene_pattern = Chem.MolFromSmarts("[C]=[C]")
    alkyne_pattern = Chem.MolFromSmarts("[C]#[C]")
    small_ring_pattern = Chem.MolFromSmarts("[r3,r4,r5,r6]")
    
    has_alcohol = mol.HasSubstructMatch(alcohol_pattern)
    has_alkene = mol.HasSubstructMatch(alkene_pattern)
    has_alkyne = mol.HasSubstructMatch(alkyne_pattern)
    has_small_ring = mol.HasSubstructMatch(small_ring_pattern)
        
    # Common VOC characteristics
    likely_voc = (
        (mol_wt < 250) or  # Small molecules
        (mol_wt < 300 and tpsa < 40) or  # Small hydrocarbons
        (mol_wt < 350 and rotatable_bonds < 12) or  # Rigid molecules
        (num_heavy_atoms < 20 and log_p < 5) or  # Small, moderately polar molecules
        (has_alcohol and mol_wt < 400) or  # Linear alcohols
        (has_alkene and mol_wt < 400) or  # Alkenes
        (has_alkyne and mol_wt < 400) or  # Alkynes
        (has_small_ring and mol_wt < 400)  # Small cyclic compounds
    )
    
    if likely_voc:
        reasons = []
        if mol_wt < 250:
            reasons.append(f"low molecular weight ({mol_wt:.1f})")
        if num_carbons < 15:
            reasons.append(f"short carbon chain (C{num_carbons})")
        if rotatable_bonds < 12:
            reasons.append("relatively rigid structure")
        if log_p < 5:
            reasons.append("moderate hydrophobicity")
        if has_alcohol:
            reasons.append("contains alcohol group")
        if has_alkene:
            reasons.append("contains alkene")
        if has_alkyne:
            reasons.append("contains alkyne")
        if has_small_ring:
            reasons.append("contains small ring")
            
        return True, "Likely VOC due to: " + ", ".join(reasons)
    
    # For borderline cases, use a more complex heuristic
    volatility_score = (
        -0.005 * mol_wt +  # Lower weight = more volatile
        -0.3 * rotatable_bonds +  # Fewer rotatable bonds = more volatile
        -0.1 * tpsa +  # Lower polarity = more volatile
        -0.2 * log_p +  # Moderate logP is preferred
        (2 if has_alcohol else 0) +  # Common VOC features
        (1 if has_alkene else 0) +
        (1 if has_alkyne else 0) +
        (1 if has_small_ring else 0)
    )
    
    if volatility_score > -12:
        return True, "Borderline VOC based on combined physicochemical properties"
        
    return False, "Combined physicochemical properties suggest non-volatile compound"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:33262',
        'name': 'volatile organic compound',
        'definition': 'Any organic compound having an initial boiling point less '
                     'than or equal to 250 degreeC (482 degreeF) measured at a '
                     'standard atmospheric pressure of 101.3 kPa.',
    }
}
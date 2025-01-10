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
    if mol_wt > 300:
        return False, f"Molecular weight ({mol_wt:.1f}) too high for VOC"
    
    # Long carbon chains typically have high boiling points
    if num_carbons > 20:
        return False, f"Carbon chain too long (C{num_carbons}) for VOC"
    
    # Many rotatable bonds indicate flexibility and higher boiling point
    if rotatable_bonds > 15:
        return False, f"Too many rotatable bonds ({rotatable_bonds}) for VOC"
    
    # High polarity (TPSA) generally means higher boiling point
    if tpsa > 100:
        return False, f"Too polar (TPSA={tpsa:.1f}) for VOC"

    # Extremely hydrophobic compounds (high logP) tend to have high boiling points
    if log_p > 6:
        return False, f"Too hydrophobic (logP={log_p:.1f}) for VOC"
        
    # Common VOC characteristics
    likely_voc = (
        (mol_wt < 200) or  # Small molecules
        (num_carbons <= 12 and tpsa < 40) or  # Small hydrocarbons
        (mol_wt < 250 and rotatable_bonds < 8) or  # Rigid small molecules
        (num_heavy_atoms < 15 and log_p < 4)  # Small, moderately polar molecules
    )
    
    if likely_voc:
        reasons = []
        if mol_wt < 200:
            reasons.append(f"low molecular weight ({mol_wt:.1f})")
        if num_carbons <= 12:
            reasons.append(f"short carbon chain (C{num_carbons})")
        if rotatable_bonds < 8:
            reasons.append("relatively rigid structure")
        if log_p < 4:
            reasons.append("moderate hydrophobicity")
            
        return True, "Likely VOC due to: " + ", ".join(reasons)
    
    # For borderline cases, use a more complex heuristic
    volatility_score = (
        -0.01 * mol_wt +  # Lower weight = more volatile
        -0.5 * rotatable_bonds +  # Fewer rotatable bonds = more volatile
        -0.2 * tpsa +  # Lower polarity = more volatile
        -0.5 * log_p  # Moderate logP is preferred
    )
    
    if volatility_score > -15:
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
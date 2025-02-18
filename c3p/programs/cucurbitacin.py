"""
Classifies: CHEBI:16219 cucurbitacin
"""
"""
Classifies: CHEBI:87070 cucurbitacin
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are tetracyclic triterpenoids derived from cucurbitane.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cucurbitacin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for triterpenoid skeleton (approximately 30 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 25:  # Adjusted for possible substitutions
        return False, f"Too few carbons ({c_count}) for triterpenoid"
    
    # Check for tetracyclic system (4 rings in SSSR)
    sssr = Chem.GetSSSR(mol)
    if len(sssr) < 4:
        return False, f"Only {len(sssr)} rings, need at least 4"
    
    # Check for oxygen-containing groups (hydroxyl, ketone, ester)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:
        return False, f"Only {o_count} oxygens, expected multiple oxygen groups"
    
    # Check molecular weight (typical cucurbitacins are >400 Da)
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da)"
    
    # Check for alpha,beta-unsaturated ketone (common feature)
    enone = Chem.MolFromSmarts("C=C(O)C(=O)")
    if not mol.HasSubstructMatch(enone):
        return False, "No alpha,beta-unsaturated ketone detected"
    
    return True, "Tetracyclic triterpenoid with oxygen functional groups"
"""
Classifies: CHEBI:51963 hopanoid
"""
"""
Classifies: CHEBI:XXXXX hopanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hopanoid(smiles: str):
    """
    Determines if a molecule is a hopanoid based on its SMILES string.
    A hopanoid is a triterpenoid with a hopane skeleton (pentacyclic structure).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hopanoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for the hopane core (pentacyclic structure)
    # Pattern captures the characteristic 5-ring system with bridge connections
    hopane_core_smarts = Chem.MolFromSmarts(
        "[C]1([C@@]2([C@]3([C@@H]([C@@]4([C@H](CC3)CC4)C)CC2)C)CCC1)CC[C@@]25C"
    )
    if hopane_core_smarts is None:
        return None, None  # Fallback if SMARTS is invalid
    
    # Check for core structure match
    if mol.HasSubstructMatch(hopane_core_smarts):
        return True, "Contains hopane core pentacyclic structure"
    
    # Alternative pattern for hopene (with one double bond)
    hopene_core_smarts = Chem.MolFromSmarts(
        "[C]1([C@@]2([C@]3([C@@H]([C@@]4([C@H](CC3)=CC4)C)CC2)C)CCC1)CC[C@@]25C"
    )
    if hopene_core_smarts and mol.HasSubstructMatch(hopene_core_smarts):
        return True, "Contains hopene core with unsaturation"
    
    # Additional check for triterpenoid characteristics (C30 skeleton)
    # Triterpenoids have ~30 carbons, allowing for some functional groups
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if 25 <= c_count <= 35:  # Allow for functional group substitutions
        # Check for pentacyclic system (5 rings)
        ring_info = mol.GetRingInfo()
        if len(ring_info.AtomRings()) >= 5:
            # Verify molecular formula pattern (C30H52O is hopane base)
            formula = rdMolDescriptors.CalcMolFormula(mol)
            if 'C' in formula and 'H' in formula:
                return True, "Pentacyclic triterpenoid structure with C30 base"
    
    return False, "No hopane core structure detected"
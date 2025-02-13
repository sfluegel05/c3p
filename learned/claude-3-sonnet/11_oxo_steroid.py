"""
Classifies: CHEBI:47787 11-oxo steroid
"""
"""
Classifies: 11-oxo steroid
Definition: Any oxo steroid that has an oxo substituent at position 11
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (True if molecule is an 11-oxo steroid, reason for classification)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for basic steroid core (four fused rings)
    # More specific pattern for steroid core with correct connectivity
    steroid_core = Chem.MolFromSmarts("[#6]1~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]3~[#6](~[#6]~[#6]~[#6]4~[#6]~[#6]~[#6]~[#6]~[#6]4~[#6]3~[#6]2)~[#6]1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # More specific SMARTS pattern for 11-oxo group in steroid context
    # The pattern ensures the ketone is at position 11 by specifying the connectivity
    # to both rings C and D of the steroid core
    oxo_11_pattern = Chem.MolFromSmarts("[#6]1~[#6]~[#6]2~[#6]~[#6]~[#6]3~[#6](=O)~[#6]~[#6]~[#6]~[#6]3~[#6]2~[#6]~[#6]1")
    
    # Find matches for the 11-oxo pattern
    matches = mol.GetSubstructMatches(oxo_11_pattern)
    if not matches:
        return False, "No ketone group at position 11"

    # Count carbons to verify it's in the typical steroid range
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 19 or carbon_count > 35:  # Expanded range to accommodate larger steroids
        return False, f"Carbon count ({carbon_count}) outside typical steroid range (19-35)"
        
    # Check molecular weight using correct descriptor module
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 600:  # Expanded range to accommodate larger steroids
        return False, f"Molecular weight ({mol_wt:.1f}) outside typical steroid range (250-600)"
        
    # Count rings
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 4:
        return False, "Insufficient number of rings for steroid structure"
    
    # Additional check for sp3 carbons typical in steroids
    sp3_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[#6](-[#6])(-[#6])(-[#6])-[#6]")))
    if sp3_count < 2:
        return False, "Insufficient sp3 carbons for typical steroid structure"

    return True, "Molecule contains steroid core with ketone group at position 11"
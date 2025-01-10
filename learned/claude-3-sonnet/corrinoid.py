"""
Classifies: CHEBI:33913 corrinoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_corrinoid(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a corrinoid based on its SMILES string.
    A corrinoid contains four reduced or partly reduced pyrrole rings joined in a 
    macrocycle by three =C- groups and one direct carbon-carbon bond linking alpha positions.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple[bool, str]: (is_corrinoid, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for basic macrocyclic framework with 4 nitrogens
    # [#7] represents nitrogen, allowing for different oxidation states
    macrocycle = Chem.MolFromSmarts('[#7]~[#6]~[#6]~[#6]~[#7]~[#6]~[#6]~[#6]~[#7]~[#6]~[#6]~[#6]~[#7]')
    if not mol.HasSubstructMatch(macrocycle):
        return False, "Missing required nitrogen macrocycle"

    # Look for four pyrrole-like rings (allowing for reduced forms)
    # This pattern matches both reduced and partially reduced pyrroles
    pyrrole = Chem.MolFromSmarts('[#7]1[#6][#6][#6][#6]1')
    matches = mol.GetSubstructMatches(pyrrole)
    if len(matches) < 4:
        return False, "Does not contain four pyrrole-like rings"

    # Look for three double bond connections between rings
    # Using (~) allows for resonance forms
    double_bond_pattern = Chem.MolFromSmarts('[#6]=[#6]')
    double_bonds = len(mol.GetSubstructMatches(double_bond_pattern))
    if double_bonds < 3:
        return False, "Missing required double bond connections"

    # Look for direct C-C bond between alpha positions
    # This pattern looks for carbons next to nitrogens that are connected
    alpha_link = Chem.MolFromSmarts('[#7][#6]-[#6][#7]')
    if not mol.HasSubstructMatch(alpha_link):
        return False, "Missing direct C-C bond between alpha positions"

    # Check ring connectivity
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings() >= 4:
        return False, "Insufficient ring count"

    # Additional check for overall size and composition
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 20:  # Corrinoids are large molecules
        return False, "Molecule too small for corrinoid structure"

    # Count nitrogens to ensure we have the core structure
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 4:
        return False, "Insufficient nitrogen atoms"

    # Success case
    reason = ("Contains four pyrrole-like rings in a macrocycle, "
             "joined by three =C- groups and one direct C-C bond")
    
    return True, reason
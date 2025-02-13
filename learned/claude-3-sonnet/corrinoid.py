"""
Classifies: CHEBI:33913 corrinoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

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

    # More specific corrin core pattern capturing the characteristic connectivity
    # Four nitrogens in specific arrangement with connecting carbons
    corrin_core = Chem.MolFromSmarts('[#7]1~[#6]~[#6]~[#6]~[#7]~[#6]~[#6]~[#6]~[#7]~[#6]~[#6]~[#6]~[#7]~[#6]~[#6]1')
    if not mol.HasSubstructMatch(corrin_core):
        return False, "Missing characteristic corrin macrocycle"

    # Check for reduced/partly reduced pyrrole rings
    # Pattern matches both fully and partially reduced forms
    pyrrole = Chem.MolFromSmarts('[#7]1(-[#6])[#6]~[#6][#6]1')
    matches = mol.GetSubstructMatches(pyrrole)
    if len(matches) < 4:
        return False, "Missing required reduced pyrrole rings"

    # Count rings to ensure macrocyclic structure
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 5:  # At least 4 pyrrole rings plus macrocycle
        return False, "Insufficient ring systems for corrinoid"

    # Check molecular composition
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if n_count < 4 or c_count < 20:
        return False, "Insufficient atoms for corrinoid structure"

    # Look for characteristic linking pattern (three =C- groups and one C-C bond)
    link_pattern = Chem.MolFromSmarts('[#7]~[#6](~[#6])~[#6]~[#7]')
    direct_link = Chem.MolFromSmarts('[#7]~[#6]-[#6]~[#7]')
    if not (mol.HasSubstructMatch(link_pattern) and mol.HasSubstructMatch(direct_link)):
        return False, "Missing characteristic linking pattern"

    # Check for cobalt (common but not required)
    has_cobalt = any(atom.GetAtomicNum() == 27 for atom in mol.GetAtoms())

    # Verify molecular weight is in typical range for corrinoids
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for corrinoid"

    # Success message varies depending on cobalt presence
    if has_cobalt:
        reason = ("Contains corrin nucleus with four modified pyrrole rings "
                 "and cobalt coordination center")
    else:
        reason = ("Contains corrin nucleus with four modified pyrrole rings "
                 "in characteristic arrangement")

    return True, reason
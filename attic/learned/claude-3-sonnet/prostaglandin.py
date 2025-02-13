"""
Classifies: CHEBI:26333 prostaglandin
"""
"""
Classifies: prostaglandin compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    Prostaglandins are naturally occurring C20 compounds derived from prostanoic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons - should be approximately 20 (allow variation for derivatives)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15 or c_count > 35:
        return False, f"Carbon count ({c_count}) outside typical range for prostaglandins"

    # More flexible cyclopentane core pattern that includes various oxidation states
    # and possible peroxide bridges (as in PGH series)
    core_patterns = [
        "[#6]1~[#6]~[#6]~[#6]~[#6]1",  # Basic 5-membered ring
        "[#6]1~[#6]~[#6](~[#8,#6])~[#6]~[#6]1",  # With substituents
        "[#6]1~[#6]2~[#6]~[#8]~[#8]~[#6]2~[#6]~[#6]~[#6]1"  # Peroxide bridge pattern
    ]
    
    ring_found = False
    for pattern in core_patterns:
        core_pat = Chem.MolFromSmarts(pattern)
        if core_pat and mol.HasSubstructMatch(core_pat):
            ring_found = True
            break
    
    if not ring_found:
        return False, "Missing characteristic cyclopentane ring structure"

    # Check for characteristic chains with double bonds
    # More flexible pattern that captures both cis and trans configurations
    chain_pattern = Chem.MolFromSmarts("[#6]~[#6]=[#6]~[#6]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Missing characteristic alkenyl chains"

    # Check for carboxylic acid or derivatives (including esters and amides)
    acid_patterns = [
        "[CX3](=[OX1])[OX2H]",  # Carboxylic acid
        "[CX3](=[OX1])[OX2][#6]",  # Ester
        "[CX3](=[OX1])[NX3]",  # Amide
        "[CX3](=[OX1])[O-]"  # Carboxylate
    ]
    
    acid_found = False
    for pattern in acid_patterns:
        acid_pat = Chem.MolFromSmarts(pattern)
        if acid_pat and mol.HasSubstructMatch(acid_pat):
            acid_found = True
            break
            
    if not acid_found:
        return False, "Missing carboxylic acid or derivative group"

    # Count oxygen atoms (allow wider range for derivatives)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, f"Insufficient oxygen atoms ({o_count}) for prostaglandin"

    # Look for oxygen substituents (more flexible pattern)
    oxygen_patterns = [
        "[#6]~[#8H1]",  # Hydroxyl
        "[#6]=[#8]",    # Ketone
        "[#8]~[#8]"     # Peroxide
    ]
    
    oxygen_found = False
    for pattern in oxygen_patterns:
        o_pat = Chem.MolFromSmarts(pattern)
        if o_pat and mol.HasSubstructMatch(o_pat):
            oxygen_found = True
            break
            
    if not oxygen_found:
        return False, "Missing characteristic oxygen substituents"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 600:  # Widened range for derivatives
        return False, f"Molecular weight ({mol_wt:.1f}) outside typical range"

    return True, "Matches prostaglandin structural features: cyclopentane ring, oxygen substituents, characteristic chains"
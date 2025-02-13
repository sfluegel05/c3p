"""
Classifies: CHEBI:26935 tetraterpenoid
"""
"""
Classifies: CHEBI:26764 tetraterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    Tetraterpenoids are derived from tetraterpenes (C40 skeleton) and may have 
    modifications like rearrangements or removal of methyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetraterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons - expanded range to include glycosylated derivatives
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 25 or c_count > 75:
        return False, f"Carbon count ({c_count}) outside typical range for tetraterpenoids (25-75)"

    # Count methyl groups - tetraterpenoids typically have multiple
    methyl_pattern = Chem.MolFromSmarts("[CH3]")
    methyl_count = len(mol.GetSubstructMatches(methyl_pattern))
    if methyl_count < 4:
        return False, f"Too few methyl groups ({methyl_count}) for tetraterpenoid"

    # Check for characteristic conjugated polyene patterns
    polyene_patterns = [
        "C=CC=CC=CC=C",  # Long conjugated chain
        "C=CC=CC=CC=CC=C",  # Extended conjugation
        "C(C)(C)=CC=C"  # Typical end group
    ]
    
    found_polyene = False
    for pattern in polyene_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_polyene = True
            break
    
    if not found_polyene:
        return False, "Missing characteristic polyene system"

    # Count double bonds - tetraterpenoids typically have many
    double_bond_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C=C")))
    if double_bond_count < 8:
        return False, f"Too few double bonds ({double_bond_count}) for tetraterpenoid"

    # Check molecular weight - adjusted range for glycosylated forms
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 1200:
        return False, f"Molecular weight ({mol_wt}) outside typical range for tetraterpenoids"

    # Count rings
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count > 6:
        return False, f"Too many rings ({ring_count}) for typical tetraterpenoid"
    
    # Count nitrogens - tetraterpenoids rarely contain nitrogen
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count > 0:
        return False, f"Contains nitrogen, unusual for tetraterpenoid"
    
    # Count oxygens - increased limit for glycosylated forms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count > 15:
        return False, f"Too many oxygens ({o_count}) for typical tetraterpenoid"

    # Check for branching - more flexible pattern
    branching_patterns = [
        "[*]([*])([*])[*]",  # General branching
        "C([C,O])([C,O])[C,O]",  # Carbon branching
        "C(=C)([C,O])[C,O]"  # Branching at double bonds
    ]
    
    found_branching = False
    for pattern in branching_patterns:
        if len(mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))) > 0:
            found_branching = True
            break
            
    if not found_branching:
        return False, "Insufficient characteristic branching"

    # Calculate degree of unsaturation
    du = rdMolDescriptors.CalcNumRotatableBonds(mol) + ring_count + double_bond_count
    if du < 8:
        return False, f"Insufficient degree of unsaturation ({du}) for tetraterpenoid"

    # Look for characteristic end groups
    end_groups = [
        "CC(C)=C",  # Typical isoprene end
        "C1C(C)=CCCC1(C)C",  # Cyclic end group
        "CC(=O)C",  # Keto end group
        "CC(O)=C"  # Hydroxy end group
    ]
    
    found_end_group = False
    for pattern in end_groups:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_end_group = True
            break
            
    if not found_end_group:
        return False, "Missing characteristic end groups"

    ring_info = f" and {ring_count} rings" if ring_count > 0 else ""
    oxygen_info = f" Contains {o_count} oxygen atoms." if o_count > 0 else ""

    return True, (f"Matches tetraterpenoid pattern with {c_count} carbons, "
                 f"{methyl_count} methyl groups, {double_bond_count} double bonds{ring_info}.{oxygen_info}")
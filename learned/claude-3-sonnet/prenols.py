"""
Classifies: CHEBI:26244 prenols
"""
"""
Classifies: prenols (CHEBI:26244)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    Prenols are alcohols with one or more isoprene units and a terminal hydroxyl or phosphate group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a prenol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for terminal groups (OH or phosphate)
    terminal_patterns = [
        Chem.MolFromSmarts("[CH2][OH]"),  # Terminal alcohol
        Chem.MolFromSmarts("[CH2]OP(=O)([O-])[O-]"),  # Terminal phosphate
        Chem.MolFromSmarts("[CH2]OP(=O)([O-])OP(=O)([O-])[O-]")  # Terminal diphosphate
    ]
    
    has_terminal = False
    terminal_type = ""
    for pattern in terminal_patterns:
        if pattern and mol.HasSubstructMatch(pattern):
            has_terminal = True
            if "P" in Chem.MolToSmiles(pattern):
                terminal_type = "phosphate"
            else:
                terminal_type = "hydroxyl"
            break
            
    if not has_terminal:
        return False, "No terminal hydroxyl or phosphate group found"

    # Count rings - prenols should be mostly linear
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count > 1:
        return False, f"Too many rings ({ring_count}) for a typical prenol"

    # Define proper isoprene unit patterns ensuring head-to-tail connectivity
    isoprene_patterns = [
        # Trans isoprene unit
        Chem.MolFromSmarts("C/C(C)=C/CC"),
        # Cis isoprene unit
        Chem.MolFromSmarts("C/C(C)=C\CC"),
        # Basic isoprene unit (any configuration)
        Chem.MolFromSmarts("CC(C)=CCC"),
        # First unit with terminal OH/phosphate
        Chem.MolFromSmarts("CC(C)=CCO"),
        Chem.MolFromSmarts("CC(C)=CCOP")
    ]
    
    # Count isoprene units
    isoprene_matches = set()
    for pattern in isoprene_patterns:
        if pattern:
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                # Add the carbons involved to avoid double counting
                isoprene_matches.update(match)
    
    num_isoprene_units = len(isoprene_matches) // 5  # Divide by typical carbons per unit
    
    if num_isoprene_units < 1:
        return False, "No complete isoprene units found"

    # Count carbons and check ratio
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    expected_c = num_isoprene_units * 5
    if abs(c_count - expected_c) > num_isoprene_units:  # Allow some deviation
        return False, f"Carbon count {c_count} not consistent with isoprene units"

    # Count double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bonds = len(mol.GetSubstructMatches(double_bond_pattern)) if double_bond_pattern else 0
    
    # Expected number of double bonds should be approximately num_isoprene_units
    if double_bonds < num_isoprene_units - 1 or double_bonds > num_isoprene_units + 1:
        return False, f"Number of double bonds ({double_bonds}) inconsistent with isoprene count"

    # Check for excessive heteroatoms
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    if n_count > 0 or s_count > 0:
        return False, "Contains unexpected heteroatoms"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Construct reason
    reason_parts = [
        f"Contains {num_isoprene_units} isoprene unit(s)",
        f"has terminal {terminal_type} group",
        f"has {double_bonds} double bond(s)",
        f"MW: {mol_wt:.1f}"
    ]
    
    return True, "; ".join(reason_parts)
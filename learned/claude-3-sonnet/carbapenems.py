"""
Classifies: CHEBI:46633 carbapenems
"""
"""
Classifies: CHEBI:23066 carbapenem
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_carbapenems(smiles: str):
    """
    Determines if a molecule is a carbapenem based on its SMILES string.
    Carbapenems have a beta-lactam ring fused to a 5-membered ring with
    substitutions at positions 2 (double bond), 3 (thioether), and 6 (hydroxyethyl).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbapenem, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core carbapenem patterns - multiple representations to catch variations
    core_patterns = [
        # Basic beta-lactam fused to 5-membered ring with C2-C3 double bond
        Chem.MolFromSmarts("[#6]1~[#6]~[#6]~[#7]2~[#6](=[O])~[#6]12"),
        # Alternative with explicit double bond
        Chem.MolFromSmarts("[#6]1-[#6]=,:[#6]-[#7]2-[#6](=O)-[#6]12"),
        # Pattern focusing on the beta-lactam fusion
        Chem.MolFromSmarts("[#6]2~[#6]~[#6]~[#7]1~[#6](=O)~[#6]1~[#6]2")
    ]
    
    core_match = False
    for pattern in core_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            core_match = True
            break
            
    if not core_match:
        return False, "Missing carbapenem core structure (fused beta-lactam and 5-membered ring)"

    # Check for beta-lactam ring (4-membered) and 5-membered ring
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    if not (4 in ring_sizes and 5 in ring_sizes):
        return False, "Must contain both 4-membered (beta-lactam) and 5-membered rings"

    # Look for key substituents and features
    # Carboxylic acid/carboxylate group at C3
    carboxyl_patterns = [
        Chem.MolFromSmarts("C(=O)O"),
        Chem.MolFromSmarts("C(=O)[O-]")
    ]
    has_carboxyl = any(mol.HasSubstructMatch(pat) for pat in carboxyl_patterns if pat is not None)

    # Common substituent patterns
    sulfur_pattern = Chem.MolFromSmarts("[#6]~[#16]")  # Carbon-Sulfur bond
    hydroxyethyl_pattern = Chem.MolFromSmarts("CC(O)")
    double_bond_pattern = Chem.MolFromSmarts("C=C")

    has_sulfur = mol.HasSubstructMatch(sulfur_pattern) if sulfur_pattern else False
    has_hydroxyethyl = mol.HasSubstructMatch(hydroxyethyl_pattern) if hydroxyethyl_pattern else False
    has_double_bond = mol.HasSubstructMatch(double_bond_pattern) if double_bond_pattern else False

    # Count key elements
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if n_count < 1:
        return False, "Must contain at least one nitrogen atom (beta-lactam)"
    if o_count < 2:
        return False, "Must contain at least two oxygen atoms"

    # Check molecular weight
    mol_wt = Chem.Descriptors.ExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for carbapenem"

    # Score features to make final decision
    features = []
    if has_carboxyl:
        features.append("carboxylic acid group")
    if has_sulfur:
        features.append("sulfur-containing substituent")
    if has_hydroxyethyl:
        features.append("hydroxyethyl group")
    if has_double_bond:
        features.append("unsaturation")

    # Require at least 2 characteristic features for positive classification
    if len(features) < 2:
        return False, "Insufficient characteristic carbapenem features"

    reason = f"Contains carbapenem core structure with {', '.join(features)}"
    return True, reason
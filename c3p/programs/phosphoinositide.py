"""
Classifies: CHEBI:18179 phosphoinositide
"""
"""
Classifies: CHEBI:18179 phosphoinositide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide is a phosphatidylinositol that is phosphorylated at one or more of the hydroxy groups of inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphoinositide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular weight - phosphoinositides are typically >600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for phosphoinositide"

    # Look for glycerol backbone pattern (C-C-C with 3 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Look for 2 ester groups (-O-C(=O)-) for fatty acid chains
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2 for fatty acid chains"

    # Look for phosphate group attached to glycerol (P-O-C)
    phosphate_glycerol_pattern = Chem.MolFromSmarts("[PX4](=[OX1])([OX2][CX4])[OX2]")
    if not mol.HasSubstructMatch(phosphate_glycerol_pattern):
        return False, "No phosphate group attached to glycerol backbone"

    # Look for inositol ring (6-membered ring with 6 oxygens, no stereochemistry requirement)
    inositol_pattern = Chem.MolFromSmarts("[C]1([C]([C]([C]([C]([C]1O)O)O)O)O)O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Look for phosphate groups attached to inositol (more flexible pattern)
    phosphate_inositol_pattern = Chem.MolFromSmarts("[PX4](=[OX1])([OX2][C]1[C]([C]([C]([C]([C]1O)O)O)O)O)[OX2]")
    phosphate_inositol_matches = mol.GetSubstructMatches(phosphate_inositol_pattern)
    
    # Must have at least one phosphate on inositol
    if len(phosphate_inositol_matches) < 1:
        return False, "No phosphorylation on inositol ring"
    
    # Count total phosphate groups (should be at least 2: one on glycerol, one on inositol)
    total_phosphates = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[PX4](=[OX1])([OX2])[OX2]")))
    if total_phosphates < 2:
        return False, f"Found only {total_phosphates} phosphate groups, need at least 2"

    # Check element counts
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    if c_count < 20:
        return False, "Too few carbons for phosphoinositide"
    if o_count < 10:
        return False, "Too few oxygens for phosphoinositide"
    if p_count < 2:
        return False, "Too few phosphorus atoms for phosphoinositide"

    return True, "Contains glycerol backbone with 2 fatty acid chains, phosphate groups, and phosphorylated inositol ring"
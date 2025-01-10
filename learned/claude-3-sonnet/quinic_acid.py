"""
Classifies: CHEBI:26493 quinic acid
"""
"""
Classifies: CHEBI:16865 quinic acid and its derivatives
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is quinic acid or a quinic acid derivative based on its SMILES string.
    Quinic acid is a cyclitol carboxylic acid with a specific substitution pattern on its cyclohexane core.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinic acid derivative, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Exclude molecules with phosphate groups (to avoid phosphatidylinositols)
    phosphate_pattern = Chem.MolFromSmarts("[P](=[O])([O,N])[O,N]")
    if mol.HasSubstructMatch(phosphate_pattern):
        return False, "Contains phosphate group - not a quinic acid"

    # Look for specific quinic acid core:
    # Cyclohexane with carboxylic acid/ester and multiple hydroxyls/esters in specific positions
    quinic_core = Chem.MolFromSmarts("""
        [C]1([C,O;!P])([$([OX2H]),$([OX2][C])])
        [C]([O,H;!P])[C]([O,H;!P])[C]([O,H;!P])[C]([O,H;!P])[C]1([O,H;!P])
    """)
    if not mol.HasSubstructMatch(quinic_core):
        return False, "No quinic acid core structure found"

    # Must have either carboxylic acid or ester
    carboxyl_pattern = Chem.MolFromSmarts("[$([CX3](=[OX1])[OX2H1]),$([CX3](=[OX1])[OX2][C])]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid or ester group found"

    # Count carbons and oxygens to ensure basic composition
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if carbon_count < 7:  # minimum carbons in quinic acid
        return False, "Too few carbons for quinic acid structure"
    
    if oxygen_count < 5:  # minimum oxygens in quinic acid
        return False, "Too few oxygens for quinic acid structure"

    # Count substitution points on the core
    core_matches = mol.GetSubstructMatches(quinic_core)
    if not core_matches:
        return False, "No valid quinic acid core found"
        
    # Get the core atoms for the first match
    core_atoms = set(core_matches[0])
    
    # Count oxygen substituents on core
    oxygen_on_core = 0
    for atom_idx in core_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:
                oxygen_on_core += 1
                
    if oxygen_on_core < 5:  # quinic acid needs at least 5 oxygens on core
        return False, "Insufficient oxygen substitution on core"

    # Look for characteristic substitution pattern
    hydroxyl_or_ester = Chem.MolFromSmarts("[$([OX2H1]),$([OX2][C](=[O])[#6])]")
    subst_count = len(mol.GetSubstructMatches(hydroxyl_or_ester))
    
    if subst_count < 4:  # Need at least 4 hydroxyl/ester groups
        return False, "Insufficient hydroxyl/ester groups"

    return True, "Contains quinic acid core structure with characteristic substitution pattern"
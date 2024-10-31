from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone (esterified prostanoid from marine corals).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a clavulone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for cyclopentenone core
    cyclopentenone_pattern = Chem.MolFromSmarts('C1CC(=O)C=C1')
    if not mol.HasSubstructMatch(cyclopentenone_pattern):
        return False, "Missing cyclopentenone core"

    # Check for ester groups
    ester_pattern = Chem.MolFromSmarts('CC(=O)O')
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "Missing acetate ester group"

    # Check for long aliphatic chain
    aliphatic_chain = Chem.MolFromSmarts('CCCCC')
    if not mol.HasSubstructMatch(aliphatic_chain):
        return False, "Missing characteristic aliphatic chain"

    # Check for at least two double bonds
    double_bonds = len(mol.GetSubstructMatches(Chem.MolFromSmarts('C=C')))
    if double_bonds < 2:
        return False, "Missing required double bonds"

    # Check for oxygen count (should have multiple oxygens due to ester groups)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 3:
        return False, "Insufficient oxygen-containing groups"

    # Optional: Check for common substituents
    halogen_pattern = Chem.MolFromSmarts('[Cl,Br,I]')
    has_halogen = mol.HasSubstructMatch(halogen_pattern)
    
    epoxide_pattern = Chem.MolFromSmarts('C1OC1')
    has_epoxide = mol.HasSubstructMatch(epoxide_pattern)

    # Check for characteristic carbon count (typically 20-25 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 15 or carbon_count > 30:
        return False, "Carbon count outside typical range for clavulones"

    # Construct reason string
    reason = "Contains cyclopentenone core, ester groups, and characteristic structure"
    if has_halogen:
        reason += " with halogen substitution"
    if has_epoxide:
        reason += " with epoxide group"

    return True, reason
# Pr=1.0
# Recall=0.9166666666666666
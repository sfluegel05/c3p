"""
Classifies: CHEBI:26267 proanthocyanidin
"""
"""
Classifies: CHEBI:26154 proanthocyanidin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule is a proanthocyanidin based on its SMILES string.
    Proanthocyanidins are flavonoid oligomers formed by condensation of hydroxyflavan units.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a proanthocyanidin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic flavonoid check - look for characteristic C6-C3-C6 skeleton
    flavonoid_pattern = Chem.MolFromSmarts("[c]1[c][c][c]([OH])[c][c]1-[CX4]-[CX4]-[c]1[c][c][c]([OH])[c][c]1")
    if not mol.HasSubstructMatch(flavonoid_pattern):
        return False, "No basic flavonoid skeleton found"

    # Look for multiple flavan units (at least 2)
    # Flavan unit pattern: hydroxyflavan core with C-ring (pyran oxygen)
    flavan_pattern = Chem.MolFromSmarts(
        "[OH]c1ccc(cc1)-[CX4]1[C@H](O)[CX4][CX4][CX4]O1")
    flavan_matches = mol.GetSubstructMatches(flavan_pattern)
    
    if len(flavan_matches) < 2:
        return False, f"Found {len(flavan_matches)} flavan units, need at least 2"

    # Check for interflavan linkages (C-C bonds between units)
    # Look for carbons connecting aromatic systems through single bonds
    linkage_pattern = Chem.MolFromSmarts("[c]-[CX4;H1,H2]-[CX4;H1,H2]-[c]")
    if not mol.HasSubstructMatch(linkage_pattern):
        return False, "No interflavan linkage pattern found"

    # Verify oligomer characteristics through molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight {mol_wt:.1f} too low for oligomer"

    # Count hydroxyl groups - proanthocyanidins typically have many -OH
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1)
    if hydroxyl_count < 4:
        return False, f"Only {hydroxyl_count} hydroxyl groups, expected ≥4"

    return True, "Contains multiple flavan units with interflavan linkages and characteristic hydroxylation"
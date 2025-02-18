"""
Classifies: CHEBI:192499 anthoxanthin
"""
"""
Classifies: Anthoxanthins (flavonoid pigments)
"""
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.
    Anthoxanthins are flavonoid pigments with a flavone core and typically glycosylation.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for flavone/flavonol core: benzopyran-4-one connected to aromatic ring
    flavone_core = Chem.MolFromSmarts("[#6]1(-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2)=,:[#6](=[#8])-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2-[#8]-1")
    if not mol.HasSubstructMatch(flavone_core):
        return False, "Missing flavone core"

    # Count hydroxyl groups (phenolic -OH)
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1)
    if hydroxyl_count < 3:
        return False, f"Only {hydroxyl_count} hydroxyl groups, need ≥3"

    # Check for glycosylation (heuristic: O-C-O or O-C-C-O patterns)
    glycoside_pattern = Chem.MolFromSmarts("[O]-[C]-[O]")  # Simple ether linkage (may miss some)
    if not mol.HasSubstructMatch(glycoside_pattern):
        # Check for common sugar patterns (e.g., hexose/pentose rings)
        sugar_ring = Chem.MolFromSmarts("[O]-[C@H]1-[C@@H](O)-[C@H](O)-[C@@H](O)-[C@H](O)-O1")
        if not mol.HasSubstructMatch(sugar_ring):
            return False, "No glycosylation detected"

    # Check oxygen count (flavones + glycosides have higher O)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 5:
        return False, f"Only {o_count} oxygens, typical anthoxanthins have ≥5"

    return True, "Flavone core with hydroxylation and glycosylation features"
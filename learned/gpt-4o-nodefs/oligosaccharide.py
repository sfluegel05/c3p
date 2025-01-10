"""
Classifies: CHEBI:50699 oligosaccharide
"""
from rdkit import Chem
from rdkit.Chem.rdchem import Atom

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    Oligosaccharides generally contain 3-10 monosaccharide units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an oligosaccharide, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of saccharide rings
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()

    # Count the number of oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Check for glycosidic bonds (C-O-C patterns)
    coc_pattern = Chem.MolFromSmarts("COC")
    coc_matches = mol.GetSubstructMatches(coc_pattern)

    # Criteria for oligosaccharide classification
    if ring_count >= 3 and 10 <= o_count <= 30 and len(coc_matches) >= 2:
        return True, "Contains features characteristic of an oligosaccharide (multiple saccharide rings and glycosidic bonds)"
    else:
        reason = (f"Characterization doesn't match oligosaccharide: "
                  f"{ring_count} rings, {o_count} oxygens, {len(coc_matches)} glycosidic-like bonds")
        return False, reason

# Example usage: Uncomment to test
# smiles_example = "OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O"
# result, reason = is_oligosaccharide(smiles_example)
# print(f"Result: {result}, Reason: {reason}")
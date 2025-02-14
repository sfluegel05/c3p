"""
Classifies: CHEBI:35381 monosaccharide
"""
"""
Classifies: CHEBI:15955 monosaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    A monosaccharide is a polyhydroxy aldehyde or polyhydroxy ketone with three or more carbon atoms,
    without glycosidic connection to other units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monosaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons, oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Must have at least 3 carbons
    if c_count < 3:
        return False, "Too few carbon atoms for monosaccharide"

    # Check for carbonyl group (aldehyde or ketone)
    carbonyl_pattern = Chem.MolFromSmarts("C(=O)")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "No carbonyl group found, not an aldehyde or ketone"

    # Count hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))

    # Check for polyhydroxy (at least 2 hydroxyls, allowing for deoxy sugars)
    if hydroxyl_count < 2:
        return False, "Not enough hydroxyl groups for monosaccharide"

    # Check for ring structures
    ring_info = mol.GetRingInfo()
    if ring_info.AtomRings():
        # Molecule contains ring structures
        # Ensure rings contain carbonyl and hydroxyl groups
        has_carbonyl_in_ring = False
        has_hydroxyl_in_ring = False
        for ring in ring_info.AtomRings():
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            ring_smarts = Chem.MolToSmarts(Chem.PathToSubmol(mol, ring_atoms, canonicalOrder=False))
            if mol.HasSubstructMatch(Chem.MolFromSmarts(f"{ring_smarts}C(=O)")):
                has_carbonyl_in_ring = True
            if mol.HasSubstructMatch(Chem.MolFromSmarts(f"{ring_smarts}[OX2H]")):
                has_hydroxyl_in_ring = True
        if not (has_carbonyl_in_ring and has_hydroxyl_in_ring):
            return False, "Ring structures do not contain carbonyl and hydroxyl groups"
    else:
        # Check for linear aldose or ketose pattern
        linear_pattern = Chem.MolFromSmarts("[CH2](C(=O))[CH](O)[CH](O)[CH](O)*")
        if not mol.HasSubstructMatch(linear_pattern):
            return False, "Does not match linear aldose or ketose pattern"

    return True, "Contains carbonyl group and at least 2 hydroxyls, meets monosaccharide criteria"
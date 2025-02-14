"""
Classifies: CHEBI:28863 flavanones
"""
"""
Classifies: Flavanones

Determines if a molecule is a flavanone based on its SMILES string.
A flavanone has a 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton.
"""

from rdkit import Chem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavanone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavanone core SMARTS pattern
    # Core structure: 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one
    flavanone_smarts = """
    [$([cR1]:[cR1]:[cR1]:[cR1]:[cR1]:[cR1])]Ar-[C@@H]1CC(=O)O[c]2[cH][cH][cH][cH]2C1
    """

    flavanone_pattern = Chem.MolFromSmarts(flavanone_smarts.strip())

    if flavanone_pattern is None:
        return None, "Error in defining flavanone SMARTS pattern"

    # Check for flavanone core
    if not mol.HasSubstructMatch(flavanone_pattern):
        # Try alternative patterns due to variability in representations
        flavanone_smarts_alt = """
        O=C1[C@H](Cc2ccccc2)Oc3ccccc13
        """
        flavanone_pattern_alt = Chem.MolFromSmarts(flavanone_smarts_alt.strip())
        if flavanone_pattern_alt is None:
            return None, "Error in defining alternative flavanone SMARTS pattern"

        if not mol.HasSubstructMatch(flavanone_pattern_alt):
            return False, "Does not contain flavanone core structure"

    # Additional checks
    # Ensure the molecule has exactly 3 rings (2 aromatic, 1 heterocyclic)
    ri = mol.GetRingInfo()
    if ri.NumRings() < 3:
        return False, "Insufficient number of rings for flavanone structure"

    aromatic_rings = 0
    heteroatom_in_ring = False
    for ring in ri.BondRings():
        ring_atoms = [mol.GetBondWithIdx(bidx).GetBeginAtom().GetIdx() for bidx in ring]
        is_aromatic = all(mol.GetAtomWithIdx(aidx).GetIsAromatic() for aidx in ring_atoms)
        if is_aromatic:
            aromatic_rings += 1
        for aidx in ring_atoms:
            if mol.GetAtomWithIdx(aidx).GetAtomicNum() not in [6, 1]:  # Not carbon or hydrogen
                heteroatom_in_ring = True
    if aromatic_rings < 2 or not heteroatom_in_ring:
        return False, "Ring structure does not match flavanone characteristics"

    return True, "Contains flavanone core structure"
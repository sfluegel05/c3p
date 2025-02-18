"""
Classifies: CHEBI:60027 polymer
"""
"""
Classifies: Polymer (CHEBI:60027)
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_polymer(smiles: str):
    """
    Determines if a molecule is a polymer based on its SMILES string.
    A polymer is a mixture of macromolecules differing in composition, length, branching, etc.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polymer, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for multiple disconnected components (mixture)
    fragments = Chem.GetMolFrags(mol, asMols=True)
    if len(fragments) > 1:
        return True, f"Mixture of {len(fragments)} components"

    # Check molecular weight (arbitrary threshold for macromolecules)
    mw = Descriptors.ExactMolWt(mol)
    if mw > 1000:
        return True, f"High molecular weight ({mw:.1f} Da)"

    # Check for repeating subunits: phosphate groups (e.g., in nucleic acids)
    phosphate = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]")
    if mol.HasSubstructMatch(phosphate):
        return True, "Contains phosphate groups indicative of polymers"

    # Check for isoprene-like repeating units (common in terpenes)
    isoprene = Chem.MolFromSmarts("C=C-C(-C)(-C)")
    isoprene_matches = len(mol.GetSubstructMatches(isoprene))
    if isoprene_matches >= 3:
        return True, f"Contains {isoprene_matches} isoprene-like units"

    # Check for multiple ester groups (e.g., polyesters)
    ester = Chem.MolFromSmarts("[#6]OC(=O)")
    ester_matches = len(mol.GetSubstructMatches(ester))
    if ester_matches >= 3:
        return True, f"Contains {ester_matches} ester groups"

    # Check for long aliphatic chains (>=10 carbons)
    long_chain = Chem.MolFromSmarts("C-C-C-C-C-C-C-C-C-C")
    if mol.HasSubstructMatch(long_chain):
        return True, "Contains a long aliphatic chain"

    # Check for polysaccharide-like patterns (glycosidic bonds)
    glycosidic = Chem.MolFromSmarts("[OX2]C1OC([C@H]([C@@H]1O)O)")
    if mol.HasSubstructMatch(glycosidic):
        return True, "Contains glycosidic bonds"

    return False, "No polymer characteristics detected"
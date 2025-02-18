"""
Classifies: CHEBI:26935 tetraterpenoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    Tetraterpenoids are derived from C40 tetraterpenes, possibly with modifications.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Basic size check - tetraterpenoids typically have ~40 carbons (allowing modifications)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30:  # Allow for some methyl removals
        return False, f"Only {c_count} carbons, insufficient for tetraterpenoid"

    # Molecular weight check (C40H64 ~536 g/mol, but modified versions may be lower)
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight {mol_wt:.1f} too low"

    # Look for isoprene (C5H8) patterns - modified to account for possible rearrangements
    # Basic isoprene unit pattern (approximate)
    isoprene_pattern = Chem.MolFromSmarts("[CH2](-[CH2])-[CH2]-[CH]=[CH2]")
    matches = len(mol.GetSubstructMatches(isoprene_pattern))
    
    # Alternative pattern for conjugated systems
    conjugated_pattern = Chem.MolFromSmarts("[CH2]=[CH]-[CH2]-[CH2]")
    conjugated_matches = len(mol.GetSubstructMatches(conjugated_pattern))
    
    # Check for at least 6 isoprene-like units (allowing for modifications)
    total_units = matches + conjugated_matches
    if total_units < 6:
        return False, f"Only {total_units} isoprene-like units found"

    # Check for methyl branches common in terpenoids
    methyl_branch = Chem.MolFromSmarts("[CH3;!R]-[CX4](-[CX4])-[CX4]")
    if len(mol.GetSubstructMatches(methyl_branch)) < 3:
        return False, "Insufficient methyl branches for terpenoid structure"

    # Check for long aliphatic chain or cyclic system with multiple double bonds
    chain_pattern = Chem.MolFromSmarts("[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No long chain structure detected"

    return True, "Meets criteria for tetraterpenoid: high carbon count, isoprene-like units, methyl branches"
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
    if c_count < 35:  # Allow for significant modifications but maintain core structure
        return False, f"Only {c_count} carbons, insufficient for tetraterpenoid"

    # Molecular weight check (C40 base ~536 g/mol, modified versions could be lower)
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 350:
        return False, f"Molecular weight {mol_wt:.1f} too low"

    # Count all double bonds - tetraterpenoids typically have multiple conjugated bonds
    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bonds < 8:
        return False, f"Only {double_bonds} double bonds, need ≥8 for conjugated system"

    # Check for extended conjugation (at least 4 consecutive double bonds)
    conjugated_system = Chem.MolFromSmarts("*=*~*=*~*=*~*=*")  # At least 4 conjugated double bonds
    if not mol.HasSubstructMatch(conjugated_system):
        return False, "No extended conjugated system detected"

    # Check for methyl groups attached to conjugated system (common in isoprenoids)
    methyl_on_conjugated = Chem.MolFromSmarts("[CH3]-[C]=[C]")
    if len(mol.GetSubstructMatches(methyl_on_conjugated)) < 3:
        return False, "Insufficient methyl groups on conjugated system"

    # Check for possible isoprene skeleton pattern (head-to-tail linkage)
    # More flexible pattern allowing modifications
    isoprene_unit = Chem.MolFromSmarts("[C]=[C]-[C]([CH3])-[CH2]-[C] |(1:5,4:0)|")
    if len(mol.GetSubstructMatches(isoprene_unit)) < 6:  # At least 6 head-tail units
        return False, "Insufficient isoprene-like linkage patterns"

    return True, "Meets criteria: high carbons, conjugated system, methyl groups, isoprene patterns"
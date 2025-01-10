"""
Classifies: CHEBI:46895 lipopeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide is a compound consisting of a peptide with attached lipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for peptide bonds (-C(=O)-N-)
    peptide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4]")
    peptide_matches = mol.GetSubstructMatches(peptide_pattern)
    if len(peptide_matches) < 2:  # Need at least 2 peptide bonds
        return False, "Insufficient peptide bonds found"

    # Look for long aliphatic chain (lipid component)
    # Match chains of at least 6 carbons
    lipid_pattern = Chem.MolFromSmarts("[CH2X4][CH2X4][CH2X4][CH2X4][CH2X4][CH2X4]")
    if not mol.HasSubstructMatch(lipid_pattern):
        return False, "No lipid chain found"

    # Count carbons and nitrogens to ensure reasonable peptide composition
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if n_count < 2:  # Need multiple nitrogens for peptide bonds
        return False, "Insufficient nitrogen atoms for peptide structure"
    
    # Calculate molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    # Typical lipopeptides have significant molecular weight and flexibility
    if mol_wt < 400:
        return False, "Molecular weight too low for typical lipopeptide"
    
    if rotatable_bonds < 10:
        return False, "Insufficient rotatable bonds for lipopeptide structure"

    # Look for amino acid characteristics (presence of NH2 or NH groups and COOH groups)
    amine_pattern = Chem.MolFromSmarts("[NX3H2,NX3H1]")
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    
    has_amine = mol.HasSubstructMatch(amine_pattern)
    has_carboxyl = mol.HasSubstructMatch(carboxyl_pattern)
    
    if not (has_amine and has_carboxyl):
        return False, "Missing characteristic amino acid groups"

    # Additional check for cyclic peptides (common in lipopeptides)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        # Look specifically for larger rings that might be peptide cycles
        large_rings = [ring for ring in ring_info.AtomRings() if len(ring) >= 6]
        if large_rings:
            return True, "Cyclic lipopeptide with attached lipid chain"

    return True, "Contains peptide bonds and lipid chain with appropriate molecular characteristics"
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

    # Look for peptide bonds (-C(=O)-N-) with connected carbons
    peptide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4]")
    peptide_matches = mol.GetSubstructMatches(peptide_pattern)
    
    if len(peptide_matches) < 3:  # Need at least 3 peptide bonds
        return False, "Insufficient peptide bonds found"

    # Look for consecutive peptide bonds (peptide chain)
    peptide_chain = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4][NX3][CX3](=[OX1])[CX4]")
    if not mol.HasSubstructMatch(peptide_chain):
        return False, "No peptide chain found"

    # Look for lipid chains - multiple patterns
    lipid_patterns = [
        "[CH2X4][CH2X4][CH2X4][CH2X4][CH2X4][CH2X4]",  # Basic chain
        "[CH2X4,CH1X4][CH2X4,CH1X4][CH2X4,CH1X4][CH2X4,CH1X4][CH2X4,CH1X4][CH2X4,CH1X4]", # Branched
        "[CH2X4,CH1X4]~[CH2X4,CH1X4]~[CH2X4,CH1X4]~[CH2X4,CH1X4]~[CH2X4,CH1X4]~[CH2X4,CH1X4]" # More general
    ]
    
    has_lipid = False
    for pattern in lipid_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_lipid = True
            break
            
    if not has_lipid:
        return False, "No lipid chain found"

    # Exclude predominantly carbohydrate structures
    sugar_pattern = Chem.MolFromSmarts("[OH1][C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) > 2:  # More than 2 sugar rings suggests glycopeptide/glycolipid
        return False, "Structure appears to be a glycopeptide or glycolipid"

    # Count key atoms and check ratios
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if n_count < 3:  # Need multiple nitrogens for peptide structure
        return False, "Insufficient nitrogen atoms for peptide structure"
    
    if o_count / n_count > 4:  # High O:N ratio suggests carbohydrate
        return False, "Atomic composition suggests glycoconjugate rather than lipopeptide"

    # Look for amino acid characteristics more specifically
    aa_patterns = [
        "[NX3H2,NX3H1][CX4H]([*])[CX3](=[OX1])[O,N]", # General amino acid
        "[NX3H2,NX3H1][CX4H](R)[CX3](=[OX1])[O,N]",   # With side chain
        "[NX3][CX4H]([*])[CX3](=[OX1])[O,N]"          # Modified amino acid
    ]
    
    aa_count = 0
    for pattern in aa_patterns:
        aa_count += len(mol.GetSubstructMatches(Chem.MolFromSmarts(pattern)))
    
    if aa_count < 2:
        return False, "Insufficient amino acid-like substructures"

    # Calculate molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    if mol_wt < 500:  # Adjusted threshold
        return False, "Molecular weight too low for typical lipopeptide"
    
    if rotatable_bonds < 15:  # Adjusted threshold
        return False, "Insufficient rotatable bonds for lipopeptide structure"

    # Check for cyclic structure
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        # Look specifically for larger rings that might be peptide cycles
        large_rings = [ring for ring in ring_info.AtomRings() if len(ring) >= 7]
        if large_rings:
            return True, "Cyclic lipopeptide with attached lipid chain"

    return True, "Contains peptide chain and lipid moiety with appropriate molecular characteristics"
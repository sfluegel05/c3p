"""
Classifies: CHEBI:25903 peptide antibiotic
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a peptide antibiotic, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check molecular weight (peptide antibiotics typically between 400-3500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for peptide antibiotic"
    if mol_wt > 3500:
        return False, "Molecular weight too high for peptide antibiotic"
    
    # Look for peptide bonds (-CO-NH-)
    amide_pattern = Chem.MolFromSmarts("[NX3H][CX3](=[OX1])")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) < 2:  # Allow smaller peptides
        return False, f"Too few peptide bonds ({len(amide_matches)}), need at least 2"
    
    # Look for specific antibiotic features
    antibiotic_features = []
    
    # Check for unusual amino acids (N-methylated, D-amino acids)
    n_methyl_pattern = Chem.MolFromSmarts("[NX3;H0]([CH3])[CX3](=[OX1])")
    if mol.HasSubstructMatch(n_methyl_pattern):
        antibiotic_features.append("N-methylated amino acids")
        
    # Look for thiazole/oxazole rings (common in peptide antibiotics)
    thiazole_pattern = Chem.MolFromSmarts("c1scnc1")
    oxazole_pattern = Chem.MolFromSmarts("c1ocnc1")
    if mol.HasSubstructMatch(thiazole_pattern) or mol.HasSubstructMatch(oxazole_pattern):
        antibiotic_features.append("contains thiazole/oxazole rings")
    
    # Look for fatty acid chains
    fatty_chain = Chem.MolFromSmarts("CCCCCCCC")  # At least 8 carbons
    if mol.HasSubstructMatch(fatty_chain):
        antibiotic_features.append("contains fatty acid chain")
    
    # Look for sugar moieties
    sugar_pattern = Chem.MolFromSmarts("[OH1][CH1]1[OH1][CH1]([CH1][CH1]1)[OH1]")
    if mol.HasSubstructMatch(sugar_pattern):
        antibiotic_features.append("contains sugar moiety")
        
    # Check for unusual elements
    has_unusual_elements = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in [16, 17]:  # S or Cl
            has_unusual_elements = True
            break
    if has_unusual_elements:
        antibiotic_features.append("contains unusual elements (S/Cl)")
    
    # Count nitrogen atoms
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 2:
        return False, "Too few nitrogen atoms for a peptide"
    
    # Calculate ratio of heteroatoms to carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    hetero_ratio = (n_count + o_count) / (c_count if c_count > 0 else 1)
    
    # Look for cyclic structure
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    
    # Calculate number of rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    # Look for basic residues
    basic_pattern = Chem.MolFromSmarts("[NX3H2,NX4H3]")
    basic_matches = mol.GetSubstructMatches(basic_pattern)
    
    # Combine evidence
    evidence = antibiotic_features.copy()
    if len(amide_matches) >= 2:
        evidence.append(f"contains {len(amide_matches)} peptide bonds")
    if ring_count > 0:
        evidence.append(f"has {ring_count} rings")
    if hetero_ratio > 0.3:
        evidence.append("high heteroatom content")
    if n_rotatable > 10:
        evidence.append("flexible structure")
    if len(basic_matches) > 0:
        evidence.append("contains basic groups")
    if 400 <= mol_wt <= 3500:
        evidence.append(f"appropriate molecular weight ({int(mol_wt)} Da)")
    
    # Make final decision
    # Must have peptide bonds AND at least one specific antibiotic feature
    if (len(amide_matches) >= 2 and             # Multiple peptide bonds
        len(antibiotic_features) >= 1 and       # At least one antibiotic feature
        hetero_ratio > 0.2):                    # Significant heteroatom content
        
        return True, f"Peptide antibiotic-like structure: {', '.join(evidence)}"
    
    return False, "Does not match peptide antibiotic characteristics"
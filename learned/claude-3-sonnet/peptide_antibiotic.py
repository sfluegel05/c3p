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
    
    # Check molecular weight (peptide antibiotics typically between 500-3500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for peptide antibiotic"
    if mol_wt > 3500:
        return False, "Molecular weight too high for peptide antibiotic"
    
    # Look for peptide bonds (-CO-NH-)
    amide_pattern = Chem.MolFromSmarts("[NX3H][CX3](=[OX1])")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) < 4:  # Require at least 4 peptide bonds
        return False, f"Too few peptide bonds ({len(amide_matches)}), need at least 4"
    
    # Initialize features list
    antibiotic_features = []
    
    # Check for cyclic peptide structure
    cyclic_peptide = Chem.MolFromSmarts("[NX3H][CX3](=[OX1])[CX4][CX4][NX3H][CX3](=[OX1])")
    if mol.HasSubstructMatch(cyclic_peptide):
        antibiotic_features.append("cyclic peptide structure")
    
    # Check for D-amino acids (look for specific stereochemistry patterns)
    d_amino = Chem.MolFromSmarts("[$([C@H]([NH2])[CX3](=O))]")
    if mol.HasSubstructMatch(d_amino):
        antibiotic_features.append("contains D-amino acids")
    
    # Check for unusual amino acids (N-methylated, etc)
    n_methyl_pattern = Chem.MolFromSmarts("[NX3;H0]([CH3])[CX3](=[OX1])")
    if mol.HasSubstructMatch(n_methyl_pattern):
        antibiotic_features.append("N-methylated amino acids")
        
    # Look for thiazole/oxazole rings
    thiazole_pattern = Chem.MolFromSmarts("c1scnc1")
    oxazole_pattern = Chem.MolFromSmarts("c1ocnc1")
    if mol.HasSubstructMatch(thiazole_pattern) or mol.HasSubstructMatch(oxazole_pattern):
        antibiotic_features.append("contains thiazole/oxazole rings")
    
    # Look for lipopeptide features (fatty acid chain + cyclic peptide)
    fatty_chain = Chem.MolFromSmarts("CCCCCCCC")  # At least 8 carbons
    if mol.HasSubstructMatch(fatty_chain) and mol.HasSubstructMatch(cyclic_peptide):
        antibiotic_features.append("lipopeptide structure")
    
    # Look for glycopeptide features but don't count as positive unless other features present
    sugar_pattern = Chem.MolFromSmarts("[OH1][CH1]1[OH1][CH1]([CH1][CH1]1)[OH1]")
    has_sugars = mol.HasSubstructMatch(sugar_pattern)
    
    # Check for unusual elements
    has_unusual_elements = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in [16, 17]:  # S or Cl
            has_unusual_elements = True
            break
    if has_unusual_elements:
        antibiotic_features.append("contains unusual elements (S/Cl)")
    
    # Count heteroatoms
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Calculate heteroatom ratio
    hetero_ratio = (n_count + o_count) / (c_count if c_count > 0 else 1)
    
    # Analyze ring structure
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    
    # Look for basic residues (common in antibiotics)
    basic_pattern = Chem.MolFromSmarts("[NX3H2,NX4H3]")
    basic_matches = mol.GetSubstructMatches(basic_pattern)
    
    # Calculate number of rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    # Combine evidence
    evidence = antibiotic_features.copy()
    if len(amide_matches) >= 4:
        evidence.append(f"contains {len(amide_matches)} peptide bonds")
    if ring_count > 0:
        evidence.append(f"has {ring_count} rings")
    if hetero_ratio > 0.3:
        evidence.append("high heteroatom content")
    if n_rotatable > 10:
        evidence.append("flexible structure")
    if len(basic_matches) > 0:
        evidence.append("contains basic groups")
    if 500 <= mol_wt <= 3500:
        evidence.append(f"appropriate molecular weight ({int(mol_wt)} Da)")
    
    # Make final decision
    # Must have multiple peptide bonds AND specific antibiotic features
    if (len(amide_matches) >= 4 and                     # Multiple peptide bonds
        (len(antibiotic_features) >= 2 or               # At least two antibiotic features
         ("cyclic peptide structure" in antibiotic_features and len(amide_matches) >= 6)) and  # Or cyclic with many peptide bonds
        hetero_ratio > 0.2 and                          # Significant heteroatom content
        not (has_sugars and len(antibiotic_features) < 2)):  # Avoid pure glycopeptides
        
        return True, f"Peptide antibiotic-like structure: {', '.join(evidence)}"
    
    return False, "Does not match peptide antibiotic characteristics"
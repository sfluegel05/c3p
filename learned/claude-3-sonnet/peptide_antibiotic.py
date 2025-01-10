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
    
    # Check molecular weight (most peptide antibiotics are >500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for peptide antibiotic"
    
    # Look for peptide bonds (-CO-NH-)
    amide_pattern = Chem.MolFromSmarts("[NX3H][CX3](=[OX1])")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) < 3:  # Most peptide antibiotics have multiple peptide bonds
        return False, f"Too few peptide bonds ({len(amide_matches)}), need at least 3"
    
    # Count nitrogen atoms (peptides typically have many N atoms)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 3:
        return False, "Too few nitrogen atoms for a peptide"
    
    # Look for amino acid-like substructures
    aa_pattern = Chem.MolFromSmarts("[NX3H][CX4H]([*])[CX3](=[OX1])[*]")
    aa_matches = mol.GetSubstructMatches(aa_pattern)
    if len(aa_matches) < 2:  # Should have multiple amino acid residues
        return False, "Insufficient amino acid-like substructures"
    
    # Check for cyclic structure (many peptide antibiotics are cyclic)
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    
    # Calculate ratio of heteroatoms to carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    hetero_ratio = (n_count + o_count) / (c_count if c_count > 0 else 1)
    
    # Calculate number of rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    # Look for basic residues (common in antimicrobial peptides)
    basic_pattern = Chem.MolFromSmarts("[NX3H2,NX4H3]")  # Primary and secondary amines
    basic_matches = mol.GetSubstructMatches(basic_pattern)
    
    # Combine evidence
    evidence = []
    if len(amide_matches) >= 3:
        evidence.append(f"contains {len(amide_matches)} peptide bonds")
    if ring_count > 0:
        evidence.append(f"has {ring_count} rings")
    if hetero_ratio > 0.3:
        evidence.append("high heteroatom content")
    if n_rotatable > 10:
        evidence.append("flexible structure")
    if len(basic_matches) > 0:
        evidence.append("contains basic groups")
    if mol_wt > 500:
        evidence.append(f"appropriate molecular weight ({int(mol_wt)} Da)")
    
    # Make final decision
    if (len(amide_matches) >= 3 and  # Multiple peptide bonds
        mol_wt > 500 and             # Substantial size
        len(aa_matches) >= 2 and     # Multiple amino acid residues
        hetero_ratio > 0.2):         # Significant heteroatom content
        
        return True, f"Peptide antibiotic-like structure: {', '.join(evidence)}"
    
    return False, "Does not match peptide antibiotic characteristics"
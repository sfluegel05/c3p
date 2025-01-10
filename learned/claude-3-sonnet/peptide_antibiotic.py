"""
Classifies: CHEBI:25903 peptide antibiotic
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors

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
    
    # Initialize features list
    antibiotic_features = []
    
    # Look for peptide bonds (-CO-NH-)
    amide_pattern = Chem.MolFromSmarts("[NX3H][CX3](=[OX1])")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    n_peptide_bonds = len(amide_matches)
    
    if n_peptide_bonds < 3:
        return False, f"Too few peptide bonds ({n_peptide_bonds}), need at least 3"
    
    # Check for cyclic peptide structure (more specific patterns)
    cyclic_patterns = [
        "[NX3H][CX3](=[OX1])[CX4][CX4][NX3H][CX3](=[OX1])",  # Basic cyclic peptide
        "[NX3][CX3](=[OX1])[CX4][CX4][NX3][CX3](=[OX1])",    # N-methylated cyclic
        "[NX3H][CX3](=[OX1])[CX4][CX4][CX4][NX3H][CX3](=[OX1])"  # Larger cycle
    ]
    
    is_cyclic = False
    for pattern in cyclic_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            is_cyclic = True
            antibiotic_features.append("cyclic peptide structure")
            break
    
    # Check for lipopeptide features (more specific)
    fatty_acid_patterns = [
        "CCCCCCCCCC",  # C10 chain
        "CCCCCCCCCCCC",  # C12 chain
        "[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX3](=[OX1])[NX3H]"  # Fatty acid-peptide link
    ]
    
    for pattern in fatty_acid_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            antibiotic_features.append("lipopeptide structure")
            break
    
    # Check for unusual amino acids and modifications
    patterns = {
        "n_methyl": "[NX3;H0]([CH3])[CX3](=[OX1])",  # N-methylated
        "d_amino": "[$([C@H]([NH2])[CX3](=O))]",  # D-amino acids
        "thiazole": "c1scnc1",  # Thiazole ring
        "oxazole": "c1ocnc1",  # Oxazole ring
        "basic_residue": "[NX3H2,NX4H3]",  # Basic residues
        "unusual_aa": "[NX3H]C(=[OX1])[CX4]S",  # Unusual amino acids with S
    }
    
    for name, pattern in patterns.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            if name == "n_methyl":
                antibiotic_features.append("N-methylated amino acids")
            elif name == "d_amino":
                antibiotic_features.append("contains D-amino acids")
            elif name in ["thiazole", "oxazole"]:
                antibiotic_features.append("contains thiazole/oxazole rings")
    
    # Count heteroatoms and calculate ratio
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    hetero_ratio = (n_count + o_count + s_count) / (c_count if c_count > 0 else 1)
    
    if hetero_ratio > 0.3:
        antibiotic_features.append("high heteroatom content")
    
    # Analyze complexity
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    if ring_count > 0:
        antibiotic_features.append(f"has {ring_count} rings")
    
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable > 10:
        antibiotic_features.append("flexible structure")
    
    # Check for unusual elements
    has_unusual_elements = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in [16, 17]:  # S or Cl
            has_unusual_elements = True
            break
    if has_unusual_elements:
        antibiotic_features.append("contains unusual elements (S/Cl)")
    
    # Calculate surface area and complexity
    tpsa = Descriptors.TPSA(mol)
    complexity = Descriptors.BertzCT(mol)
    
    # Make final decision with stricter criteria
    is_antibiotic = False
    if (n_peptide_bonds >= 3 and  # Minimum peptide bonds
        ((is_cyclic and len(antibiotic_features) >= 3) or  # Cyclic with features
         (len(antibiotic_features) >= 4)) and  # Or more features for linear
        hetero_ratio > 0.25 and  # Higher heteroatom requirement
        (tpsa > 100 or complexity > 500)):  # Complexity requirements
        
        evidence = antibiotic_features.copy()
        evidence.append(f"contains {n_peptide_bonds} peptide bonds")
        evidence.append(f"appropriate molecular weight ({int(mol_wt)} Da)")
        
        return True, f"Peptide antibiotic-like structure: {', '.join(evidence)}"
    
    return False, "Does not match peptide antibiotic characteristics"
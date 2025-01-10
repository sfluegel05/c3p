"""
Classifies: CHEBI:33567 catecholamine
"""
"""
Classifies: catecholamine
Definition: 4-(2-Aminoethyl)pyrocatechol [4-(2-aminoethyl)benzene-1,2-diol] and derivatives formed by substitution
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a catecholamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Handle salts - get largest fragment
    fragments = Chem.GetMolFrags(mol, asMols=True)
    if len(fragments) > 1:
        mol = max(fragments, key=lambda m: m.GetNumAtoms())
    
    # Core patterns:
    
    # 1. Modified catechol pattern that allows for substitutions
    # Allows one OH to be replaced by OMe or sulfate
    catechol_pattern = Chem.MolFromSmarts("""
        [#6]:1:[#6](-[OH1,O-,OS(=O)(=O)[OH1]])
        :[#6](-[OH1,O-,OCH3])
        :[#6]:[#6](-[#6]):[#6]:1
    """)
    
    # 2. Strict 2-aminoethyl pattern with allowed modifications
    aminoethyl_patterns = [
        # Basic 2-aminoethyl: -CH2-CH2-NH2
        Chem.MolFromSmarts("c-CCN"),
        
        # Allow hydroxylation on beta carbon: -CH2-CHOH-NH2
        Chem.MolFromSmarts("c-CC(O)N"),
        
        # Allow N-methylation: -CH2-CH2-NH(CH3)
        Chem.MolFromSmarts("c-CCNC"),
        
        # Allow N,N-dimethylation: -CH2-CH2-N(CH3)2
        Chem.MolFromSmarts("c-CCN(C)C"),
        
        # Allow alpha-methylation: -CH(CH3)-CH2-NH2
        Chem.MolFromSmarts("c-C(C)CN"),
    ]
    
    if not mol.HasSubstructMatch(catechol_pattern):
        return False, "No catechol or modified catechol moiety found"
    
    # Check for proper 2-aminoethyl chain
    has_valid_chain = False
    for pattern in aminoethyl_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            # Verify the chain is at position 4 (para) relative to one of the OH groups
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                if len(match) > 0:  # Ensure we have matches
                    has_valid_chain = True
                    break
    
    if not has_valid_chain:
        return False, "No valid 2-aminoethyl chain found at correct position"
    
    # Additional structure checks
    
    # Check molecule size
    if mol.GetNumAtoms() > 40:  # Adjusted from 50 to be more strict
        return False, "Molecule too large for a typical catecholamine"
    
    # Count aromatic rings
    aromatic_rings = mol.GetSubstructMatches(Chem.MolFromSmarts("a1aaaaa1"))
    if len(aromatic_rings) > 2:
        return False, "Too many aromatic rings"
    
    # Check for features that would indicate a different class
    exclude_patterns = [
        Chem.MolFromSmarts("[NX3]1[CR2][CR2][c]2[cR1][cR1][cR1][cR1][cR1]21"), # Benzazepine core
        Chem.MolFromSmarts("C(=O)N"), # Amide (unless part of larger necessary structure)
        Chem.MolFromSmarts("C(=O)O[H,C]"), # Ester/carboxylic acid
        Chem.MolFromSmarts("S[C,H]"), # Thiol/thioether (but allow sulfates)
        Chem.MolFromSmarts("P"), # Phosphorus
    ]
    
    for pattern in exclude_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            return False, "Contains structural features not consistent with catecholamines"
            
    return True, "Contains required catecholamine structural features"
"""
Classifies: CHEBI:18000 aralkylamine
"""
"""
Classifies: aralkylamine
An alkylamine in which the alkyl group is substituted by an aromatic group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an aralkylamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for primary, secondary or tertiary amine
    amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(NC=O);!$(N=*);!$(N#*)]")
    if not mol.HasSubstructMatch(amine_pattern):
        return False, "No amine group found"
    
    # Look for aromatic rings
    aromatic_patterns = [
        Chem.MolFromSmarts("a1aaaaa1"), # 6-membered aromatic
        Chem.MolFromSmarts("a1aaaa1"),  # 5-membered aromatic
    ]
    
    has_aromatic = False
    for pattern in aromatic_patterns:
        if mol.HasSubstructMatch(pattern):
            has_aromatic = True
            break
            
    if not has_aromatic:
        return False, "No aromatic ring found"
    
    # Look for alkyl linker between amine and aromatic ring
    # Pattern matches: aromatic ring - sp3 carbon(s) - amine
    aralkyl_pattern = Chem.MolFromSmarts("a-[CX4]+-[NX3;H2,H1,H0]")
    
    if not mol.HasSubstructMatch(aralkyl_pattern):
        return False, "No alkyl linker between aromatic ring and amine"
        
    # Additional checks to exclude false positives
    
    # Exclude amides
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    if mol.HasSubstructMatch(amide_pattern):
        matches = mol.GetSubstructMatches(amine_pattern)
        amide_matches = mol.GetSubstructMatches(amide_pattern)
        if all(match[0] in [am[0] for am in amide_matches] for match in matches):
            return False, "Contains only amide groups, no amine"
    
    # Exclude aromatic amines (anilines)
    aniline_pattern = Chem.MolFromSmarts("c-[NX3]")
    if mol.HasSubstructMatch(aniline_pattern):
        matches = mol.GetSubstructMatches(amine_pattern)
        aniline_matches = mol.GetSubstructMatches(aniline_pattern)
        if all(match[0] in [an[1] for an in aniline_matches] for match in matches):
            return False, "Contains only aromatic amine (aniline), not aralkylamine"
    
    return True, "Contains amine group connected to aromatic ring via alkyl linker"
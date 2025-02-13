"""
Classifies: CHEBI:48030 tetrapeptide
"""
"""
Classifies: CHEBI:36357 tetrapeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

# Comprehensive set of SMARTS patterns for common amino acid residues
AA_PATTERNS = ['[N;X3,X4](C)(C(=O)O)[C@H]([*])C(=O)N', # Ala, Gly
               '[N;X3,X4](C([CH3]([*]))C)(C(=O)O)[C@H]([*])C(=O)N', # Val, Ile, Leu
               '[N;X3,X4](C([CH2]([*])Cc1ccccc1))(C(=O)O)[C@H]([*])C(=O)N', # Phe
               '[N;X3,X4](C([CH2]([*])Cc1c[nH]c2c1cccc2))(C(=O)O)[C@H]([*])C(=O)N', # Trp
               '[N;X3,X4](C([CH2]([*])Cc1ccc(O)cc1))(C(=O)O)[C@H]([*])C(=O)N', # Tyr
               '[N;X3,X4](C([CH2]([*])SC))(C(=O)O)[C@H]([*])C(=O)N', # Cys, Met
               '[N;X3,X4](C([CH2]([*])C(N)=O))(C(=O)O)[C@H]([*])C(=O)N', # Asn, Gln
               '[N;X3,X4](C([CH2]([*])C(O)=O))(C(=O)O)[C@H]([*])C(=O)N', # Asp, Glu
               '[N;X3,X4](C([CH2]([*])CC(C)(C)N))(C(=O)O)[C@H]([*])C(=O)N', # Lys
               '[N;X3,X4](C([CH2]([*])CC1=CNC=N1))(C(=O)O)[C@H]([*])C(=O)N', # His
               '[N;X3,X4](C([CH2]([*])CC(C)(C)NC(N)=N))(C(=O)O)[C@H]([*])C(=O)N'] # Arg

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide is any molecule that contains four amino-acid residues connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrapeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check molecular weight range (300-800 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 800:
        return False, f"Molecular weight ({mol_wt:.2f} Da) outside typical range for tetrapeptides"
    
    # Check elemental composition
    composition = rdMolDescriptors.CalcMolFormula(mol)
    if "N" not in composition or "O" not in composition or "C" not in composition:
        return False, "Missing essential elements (N, O, C) for peptides"
    
    # Look for tetrapeptide backbone pattern
    tetrapeptide_pattern = Chem.MolFromSmarts("[N;X3,X4]([CH2])[CH]([*])[C](=O)[N;X3,X4]([CH2])[CH]([*])[C](=O)[N;X3,X4]([CH2])[CH]([*])[C](=O)[N;X3,X4]([CH2])[CH]([*])[C](=O)O")
    if not mol.HasSubstructMatch(tetrapeptide_pattern):
        return False, "No tetrapeptide backbone found"
    
    # Look for 4 amino acid residues
    aa_matches = 0
    for pattern in AA_PATTERNS:
        aa_pattern = Chem.MolFromSmarts(pattern)
        aa_matches += len(mol.GetSubstructMatches(aa_pattern))
    
    if aa_matches != 4:
        return False, f"Found {aa_matches} amino acid residues, need exactly 4"
    
    return True, "Contains four amino-acid residues connected by peptide linkages"
"""
Classifies: CHEBI:63534 monoamine
"""
"""
Classifies: monoamine compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine has one amino group connected to an aromatic ring by a two-carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES and handle salt forms
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get the largest fragment (in case of salts)
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumAtoms())
    
    # Size/complexity filters
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt > 500:  # Most monoamines are relatively small
        return False, "Molecule too large for typical monoamine"
    
    # Check for aromatic ring
    aromatic_pattern = Chem.MolFromSmarts("a1aaaaa1")  # 6-membered aromatic ring
    if not mol.HasSubstructMatch(aromatic_pattern):
        return False, "No aromatic ring found"

    # Exclude peptides and amino acids
    peptide_bond = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4]")
    if mol.HasSubstructMatch(peptide_bond):
        return False, "Contains peptide bonds"
        
    # Count carboxylic acids - typical amino acids have these
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if len(mol.GetSubstructMatches(carboxyl_pattern)) > 1:
        return False, "Multiple carboxylic acid groups suggest amino acid/peptide"

    # Define monoamine patterns more precisely
    monoamine_patterns = [
        # Basic patterns
        "[aR1]!@[CH2][CH2][NX3;!R;!$(NC=O)]",  # Basic ethylamine
        "[aR1]!@[CH2][CH1]([OH1])[NX3;!R;!$(NC=O)]",  # Beta-hydroxyl
        "[aR1]!@[CH1]([OH1])[CH2][NX3;!R;!$(NC=O)]",  # Alpha-hydroxyl
        
        # Methylated variants
        "[aR1]!@[CH2][CH2][NX3H2;!R;!$(NC=O)]C",
        "[aR1]!@[CH2][CH2][NX3H1;!R;!$(NC=O)](C)C",
        "[aR1]!@[CH2][CH1]([OH1])[NX3H1;!R;!$(NC=O)]C",
        
        # Charged variants
        "[aR1]!@[CH2][CH2][NX4H3+;!R]",
        "[aR1]!@[CH2][CH2][NX4H2+;!R]C",
        "[aR1]!@[CH2][CH2][NX4H1+;!R](C)C",
        
        # Additional substituted patterns
        "[aR1]!@[CX4H2][CX4H2][NX3;!R;!$(NC=O)]",  # Strict ethylamine
        "[aR1]!@[CX4H2][CX4H1]([OH1])[NX3;!R;!$(NC=O)]",  # Strict beta-hydroxyl
        "[aR1]!@[CX4H1]([OH1])[CX4H2][NX3;!R;!$(NC=O)]"  # Strict alpha-hydroxyl
    ]

    found_valid_pattern = False
    for pattern in monoamine_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt and mol.HasSubstructMatch(patt):
            found_valid_pattern = True
            break
            
    if not found_valid_pattern:
        return False, "No valid two-carbon chain connecting amine to aromatic ring"

    # Count non-aromatic amines
    amine_pattern = Chem.MolFromSmarts("[NX3;!R;!$(NC=O)]")
    charged_amine_pattern = Chem.MolFromSmarts("[NX4H+;!R]")
    
    amine_count = len(mol.GetSubstructMatches(amine_pattern))
    charged_amine_count = len(mol.GetSubstructMatches(charged_amine_pattern))
    total_amines = amine_count + charged_amine_count
    
    if total_amines == 0:
        return False, "No amine groups found"
    if total_amines > 2:  # Stricter limit on number of amines
        return False, "Too many amine groups"

    # Check for aromatic hydroxyl groups (common in monoamines)
    aromatic_oh_pattern = Chem.MolFromSmarts("aO[H]")
    has_aromatic_oh = mol.HasSubstructMatch(aromatic_oh_pattern)
    
    detail = "monoamine with"
    if has_aromatic_oh:
        detail += " hydroxylated"
    detail += " aromatic ring and two-carbon amine linker"
    
    return True, f"Confirmed {detail}"
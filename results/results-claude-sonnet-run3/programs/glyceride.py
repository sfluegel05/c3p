from rdkit import Chem
from rdkit.Chem import AllChem

def is_glyceride(smiles: str):
    """
    Determines if a molecule is a glyceride (ester of glycerol with fatty acids).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a glyceride, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for glycerol backbone (3 carbons connected with oxygens)
    glycerol_pattern = Chem.MolFromSmarts("[C;!$(C=O)][C;!$(C=O)][C;!$(C=O)]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Get glycerol backbone atoms
    glycerol_match = mol.GetSubstructMatch(glycerol_pattern)

    # Check that each carbon in glycerol backbone has an oxygen attached
    o_count = 0
    for c_idx in glycerol_match:
        atom = mol.GetAtomWithIdx(c_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == 'O':
                o_count += 1

    if o_count < 3:
        return False, "Glycerol backbone missing oxygen substituents"

    # Check for ester groups or ester-like patterns with wildcards
    ester_pattern = Chem.MolFromSmarts("[#6,*]C(=O)O[#6]")
    wildcard_pattern = Chem.MolFromSmarts("[*]O[CH2,CH]O[*]")
    
    if mol.HasSubstructMatch(ester_pattern):
        ester_count = len(mol.GetSubstructMatches(ester_pattern))
        if ester_count >= 3:
            return True, "Triglyceride"
        elif ester_count == 2:
            return True, "Diglyceride"
        else:
            return True, "Monoglyceride"
    elif mol.HasSubstructMatch(wildcard_pattern):
        # Count wildcard substituents
        wildcard_count = len(mol.GetSubstructMatches(wildcard_pattern))
        if wildcard_count >= 3:
            return True, "Triglyceride with wildcard substituents"
        elif wildcard_count == 2:
            return True, "Diglyceride with wildcard substituents"
        else:
            return True, "Monoglyceride with wildcard substituents"
            
    # Check for simple glyceride pattern with wildcards
    simple_pattern = Chem.MolFromSmarts("[O,*]CC(O[*])CO[*]")
    if mol.HasSubstructMatch(simple_pattern):
        return True, "Glyceride with wildcard substituents"

    return False, "No glyceride pattern found"
# Pr=1.0
# Recall=0.8333333333333334
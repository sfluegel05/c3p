"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
"""
Classifies: nucleoside phosphate
A nucleobase-containing molecular entity that is a nucleoside in which one or more 
of the sugar hydroxy groups has been converted into a mono- or poly-phosphate.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a nucleoside phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for phosphate group
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"
    
    # Check for nucleobase patterns
    # Purine bases (adenine, guanine)
    purine_pattern = Chem.MolFromSmarts("c1ncnc2[nX3]cnc12")
    # Pyrimidine bases (cytosine, thymine, uracil)
    pyrimidine_pattern = Chem.MolFromSmarts("[nX3]1c([nX3H])cc(=O)[nX3]c1=O")
    
    has_purine = mol.HasSubstructMatch(purine_pattern)
    has_pyrimidine = mol.HasSubstructMatch(pyrimidine_pattern)
    
    if not (has_purine or has_pyrimidine):
        return False, "No nucleobase (purine or pyrimidine) found"
    
    # Check for sugar (ribose/deoxyribose) pattern
    # This pattern looks for the furanose ring connected to a nucleobase
    sugar_pattern = Chem.MolFromSmarts("[CH2X4]-[CH1X4]-[CH1X4]-[CH1X4]-[OX2]")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar moiety found"
    
    # Check connectivity: phosphate should be attached to sugar
    sugar_phosphate_pattern = Chem.MolFromSmarts("[CH2X4]-[CH1X4]-[CH1X4]-[CH1X4]-[OX2]-P(=O)([OX2])[OX2]")
    if not mol.HasSubstructMatch(sugar_phosphate_pattern):
        # Try alternate pattern for 2' or 3' phosphates
        alt_sugar_phosphate = Chem.MolFromSmarts("[CH1X4](-[OX2]-P(=O)([OX2])[OX2])-[CH1X4]-[CH1X4]-[OX2]")
        if not mol.HasSubstructMatch(alt_sugar_phosphate):
            return False, "Phosphate group not properly connected to sugar"
    
    # Count phosphates
    phosphate_matches = len(mol.GetSubstructMatches(phosphate_pattern))
    base_type = "purine" if has_purine else "pyrimidine"
    
    return True, f"Contains {base_type} nucleobase, sugar moiety, and {phosphate_matches} phosphate group(s)"
"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: CHEBI:36976 nucleotide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    A nucleotide must have a nucleoside base, sugar moiety, and at least one phosphate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleotide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphate group
    phosphate_pattern = Chem.MolFromSmarts("[PX4](=[OX1])([$([OX2H]),$([OX2-])])[$([OX2H]),$([OX2-])][$([OX2H]),$([OX2-]),$(O)]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for sugar (furanose) ring
    furanose_pattern = Chem.MolFromSmarts("[CH2]-[CH1]-[CH1]-[CH1]-O1")
    deoxy_furanose_pattern = Chem.MolFromSmarts("[CH2]-[CH2]-[CH1]-[CH1]-O1")
    has_furanose = mol.HasSubstructMatch(furanose_pattern)
    has_deoxy_furanose = mol.HasSubstructMatch(deoxy_furanose_pattern)
    
    if not (has_furanose or has_deoxy_furanose):
        return False, "No sugar (ribose/deoxyribose) moiety found"

    # Check for nucleobase patterns
    # Purine bases (adenine, guanine, or derivatives)
    purine_pattern = Chem.MolFromSmarts("c12ncnc1ncn2")
    
    # Pyrimidine bases (cytosine, thymine, uracil, or derivatives)
    pyrimidine_pattern = Chem.MolFromSmarts("c1cn([*])c(=O)nc1")
    
    has_purine = mol.HasSubstructMatch(purine_pattern)
    has_pyrimidine = mol.HasSubstructMatch(pyrimidine_pattern)
    
    if not (has_purine or has_pyrimidine):
        return False, "No nucleobase (purine/pyrimidine) found"

    # Check connection between sugar and phosphate
    # This is a simplified check - looking for C-O-P pattern
    sugar_phosphate_pattern = Chem.MolFromSmarts("[CH2]OP(=O)")
    if not mol.HasSubstructMatch(sugar_phosphate_pattern):
        return False, "Phosphate not properly connected to sugar"

    # Check connection between sugar and base
    # Looking for N-glycosidic bond
    glycosidic_pattern = Chem.MolFromSmarts("[NX3]1-[#6]-[#6]-[#6]-[#6](-[#6])-O1")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No proper N-glycosidic bond between base and sugar"

    return True, "Contains nucleobase, sugar, and phosphate with correct connectivity"
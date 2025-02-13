"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
"""
Classifies: CHEBI:63857 nucleoside phosphate

A nucleoside phosphate is a nucleobase-containing molecular entity that is a nucleoside 
in which one or more of the sugar hydroxy groups has been converted into a mono- or poly-phosphate.
The term includes both nucleotides and non-nucleotide nucleoside phosphates.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdqueries

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
    
    # Check for nucleobase
    nucleobase_pattern = Chem.MolFromSmarts("a1aa[n;r5,r6]c2nc[n;r5,r6]c(N)nc12")
    has_nucleobase = mol.HasSubstructMatch(nucleobase_pattern)
    
    # Check for sugar ring
    sugar_pattern = Chem.MolFromSmarts("OC1C(O)C(O)C(O)C1")
    has_sugar = mol.HasSubstructMatch(sugar_pattern)
    
    # Check for phosphate group(s)
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(O)=O")
    has_phosphate = mol.HasSubstructMatch(phosphate_pattern)
    
    # Check for glycosidic bond
    glycosidic_pattern = Chem.MolFromSmarts("Cn1cnc2c1ncnc2OC")
    has_glycosidic = mol.HasSubstructMatch(glycosidic_pattern)
    
    if has_nucleobase and has_sugar and has_phosphate and has_glycosidic:
        return True, "Contains nucleobase, sugar ring, phosphate group(s), and glycosidic bond"
    else:
        missing = []
        if not has_nucleobase:
            missing.append("nucleobase")
        if not has_sugar:
            missing.append("sugar ring")
        if not has_phosphate:
            missing.append("phosphate group(s)")
        if not has_glycosidic:
            missing.append("glycosidic bond")
        return False, f"Missing {', '.join(missing)} for nucleoside phosphate"
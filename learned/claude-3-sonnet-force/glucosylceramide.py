"""
Classifies: CHEBI:36500 glucosylceramide
"""
"""
Classifies: CHEBI:18099 glucosylceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    A glucosylceramide is a ceramide with a glucose moiety attached via a glycosidic or ester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosylceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for ceramide backbone patterns
    ceramide_patterns = [
        Chem.MolFromSmarts("[NX3H2;$(NC(=O)C)][CX4H](C[CX4])([CX4])[CX3](=O)[OX2H]"),  # Standard ceramide
        Chem.MolFromSmarts("[NX3H2;$(NC(=O)C)][CX4H](C[CX3]=[CX3])([CX4])[CX3](=O)[OX2H]"),  # Ceramide with double bond
        Chem.MolFromSmarts("[NX3H2;$(NC(=O)C)][CX4H](C[CX4])([CX4])[CX3](=O)[OX2H,OX1H0-]")  # Ceramide with substituents
    ]
    ceramide_match = any(mol.HasSubstructMatch(pattern) for pattern in ceramide_patterns)
    if not ceramide_match:
        return False, "No ceramide backbone found"
    
    # Look for glucose moiety patterns
    glucose_patterns = [
        Chem.MolFromSmarts("[OX2;$([C@H]1[C@@H]([C@@H]([C@@H]([C@H](O1)O)O)O)O[CX4;!$(NC=O)])][CX4]"),  # Glycosidic bond
        Chem.MolFromSmarts("[OX2;$([C@H]1[C@@H]([C@@H]([C@@H]([C@H](O1)O)O)O)O[CX3]=O)]")  # Ester bond
    ]
    glucose_match = any(mol.HasSubstructMatch(pattern) for pattern in glucose_patterns)
    if not glucose_match:
        return False, "No glucose moiety found"
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 18:
        return False, "Too few carbons for glucosylceramide"
    if o_count < 7:
        return False, "Too few oxygens for glucosylceramide"
    
    # Check molecular weight - glucosylceramides typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for glucosylceramide"
    
    return True, "Contains a ceramide backbone with a glucose moiety attached via a glycosidic or ester bond"
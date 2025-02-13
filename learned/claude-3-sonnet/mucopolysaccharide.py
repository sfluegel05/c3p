"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
"""
Classifies: CHEBI:26611 mucopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide based on its SMILES string.
    A mucopolysaccharide is a polysaccharide composed of alternating uronic acids and glycosamines,
    and commonly partially esterified with sulfuric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mucopolysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of uronic acid and glycosamine units
    has_uronic_acid = any(atom.GetHybridizationStatePerBondRingCt() == Chem.HybridizationType.SP2
                          and atom.GetTotalNumHs() == 1 for atom in mol.GetAtoms())
    has_glycosamine = any(atom.GetAtomicNum() == 7 and atom.GetTotalNumHs() == 2 for atom in mol.GetAtoms())
    
    if not has_uronic_acid or not has_glycosamine:
        return False, "Missing uronic acid or glycosamine units"
    
    # Look for sulfate groups or sulfur atoms
    has_sulfate = any(atom.GetAtomicNum() == 16 for atom in mol.GetAtoms())
    
    # Check molecular weight and atom counts
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for mucopolysaccharide"
    
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 10 or o_count < 5:
        return False, "Insufficient carbon and oxygen atoms for polysaccharide"
    
    # Classify as mucopolysaccharide if uronic acid, glycosamine, and sulfation present
    if has_uronic_acid and has_glycosamine and has_sulfate:
        return True, "Contains uronic acid and glycosamine units, partially esterified with sulfate groups"
    else:
        return False, "Does not match mucopolysaccharide structural features"
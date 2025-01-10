"""
Classifies: CHEBI:36092 clavulone
"""
"""
Classifies: clavulone
Definition: A class of esterified prostanoids obtained from marine corals.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a clavulone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for cyclopentenone core
    cyclopentenone = Chem.MolFromSmarts("[#6]1[#6][#6][#6](=[O])[#6]1")
    if not mol.HasSubstructMatch(cyclopentenone):
        return False, "No cyclopentenone core found"
    
    # Check for at least one ester group
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[#6]")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group found"
    
    # Check for conjugated double bonds
    conjugated_db = Chem.MolFromSmarts("C=CC=C")
    if not mol.HasSubstructMatch(conjugated_db):
        return False, "No conjugated double bonds found"
    
    # Check for long carbon chain (at least 5 carbons in sequence)
    long_chain = Chem.MolFromSmarts("[#6]~[#6]~[#6]~[#6]~[#6]")
    if not mol.HasSubstructMatch(long_chain):
        return False, "No long carbon chain found"
    
    # Count carbons and oxygens to ensure reasonable size
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 15:
        return False, "Too few carbons for clavulone"
    if o_count < 2:
        return False, "Too few oxygens for clavulone"
    
    # Check molecular weight - clavulones typically >300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for clavulone"
    
    # Check for specific substitution pattern around cyclopentenone
    # Looking for a carbon chain attached to the ring
    specific_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#6](=[O])[#6]1[#6,O]")
    if not mol.HasSubstructMatch(specific_pattern):
        return False, "Incorrect substitution pattern around cyclopentenone"
    
    # Count number of double bonds
    double_bonds = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C=C")))
    if double_bonds < 2:
        return False, "Too few double bonds for clavulone"
    
    return True, "Contains cyclopentenone core with characteristic substitution pattern, ester groups, and conjugated system"
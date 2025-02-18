"""
Classifies: CHEBI:27325 xanthophyll
"""
"""
Classifies: CHEBI:27345 xanthophyll
"""
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    Xanthophylls are oxygenated carotenoids with a conjugated system and oxygen-containing groups.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Verify carotenoid-like structure with extended conjugated system
    polyene_pattern = Chem.MolFromSmarts("C=CC=C.C=CC=C.C=CC=C")  # At least 3 conjugated double bond pairs
    if not mol.HasSubstructMatch(polyene_pattern):
        return False, "No carotenoid backbone"
    
    # Check for oxygen atoms
    if sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 8) < 1:
        return False, "No oxygen"
    
    # Detect oxygen-containing groups (improved patterns)
    hydroxyl = mol.HasSubstructMatch(Chem.MolFromSmarts("[OH]"))
    epoxy = mol.HasSubstructMatch(Chem.MolFromSmarts("[C;!$(C=O)]O[C;!$(C=O)]"))  # Better epoxy match
    ketone = mol.HasSubstructMatch(Chem.MolFromSmarts("C=O"))
    ester = mol.HasSubstructMatch(Chem.MolFromSmarts("[#6]OC(=O)"))
    ether = mol.HasSubstructMatch(Chem.MolFromSmarts("[CX4]O[CX4]"))
    
    if not any([hydroxyl, epoxy, ketone, ester, ether]):
        return False, "No oxygen functional groups"
    
    # Check conjugated system length in main chain
    conjugated_system = 0
    max_conj = 0
    for bond in mol.GetBonds():
        if bond.GetIsConjugated():
            conjugated_system += 1
        else:
            max_conj = max(max_conj, conjugated_system)
            conjugated_system = 0
    max_conj = max(max_conj, conjugated_system)
    
    if max_conj < 6:  # Allow shorter systems with oxygen groups
        return False, f"Conjugated system too short ({max_conj} bonds)"
    
    # Molecular weight check (relaxed)
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da)"
    
    return True, "Oxygenated carotenoid with conjugated system"
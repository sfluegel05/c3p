"""
Classifies: CHEBI:36249 bile acid conjugate
"""
"""
Classifies: CHEBI bile acid conjugate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    Bile acid conjugates have a steroid core conjugated to groups like glycine, taurine, sulfate, etc.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Define conjugated group patterns
    glycine = Chem.MolFromSmarts("[NH1]C(=O)CC(=O)O")  # Glycine conjugate
    taurine = Chem.MolFromSmarts("[NH1]C(=O)CCS(=O)(=O)O")  # Taurine conjugate
    sulfate = Chem.MolFromSmarts("[OX2]S(=O)(=O)[OX1]")  # Sulfate ester
    glucuronate = Chem.MolFromSmarts("[OX2]C(=O)C1C(O)C(O)C(O)C(O1)")  # Glucuronate ester (simplified)
    
    # Check for any conjugated groups
    conjugated = []
    if mol.HasSubstructMatch(glycine):
        conjugated.append("glycine")
    if mol.HasSubstructMatch(taurine):
        conjugated.append("taurine")
    if mol.HasSubstructMatch(sulfate):
        conjugated.append("sulfate")
    if mol.HasSubstructMatch(glucuronate):
        conjugated.append("glucuronate")
    
    if not conjugated:
        return False, "No conjugated groups detected"
    
    # Check for steroid-like features (simplified)
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 3:
        return False, f"Insufficient rings ({n_rings}) for steroid core"
    
    hydroxyls = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() > 0)
    if hydroxyls < 1:
        return False, "No hydroxyl groups found"
    
    # Check molecular weight (typical bile acid conjugates are >400 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da)"
    
    return True, f"Conjugated groups: {', '.join(conjugated)} with steroid features"
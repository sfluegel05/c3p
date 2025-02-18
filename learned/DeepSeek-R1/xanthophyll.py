"""
Classifies: CHEBI:27325 xanthophyll
"""
"""
Classifies: CHEBI:27345 xanthophyll
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    Xanthophylls are oxygenated carotenoids with a long conjugated system and oxygen-containing groups.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check oxygen presence
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count == 0:
        return False, "No oxygen atoms"
    
    # Check for minimum number of conjugated double bonds (carotenoid characteristic)
    conjugated_double = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetIsConjugated():
            conjugated_double += 1
    if conjugated_double < 8:
        return False, f"Only {conjugated_double} conjugated double bonds"
    
    # Check for oxygen-containing functional groups (hydroxyl, epoxy, ketone, ester)
    hydroxyl = Chem.MolFromSmarts("[OH]")
    epoxy = Chem.MolFromSmarts("C1OC1")  # Simple epoxy pattern
    ketone = Chem.MolFromSmarts("C=O")
    ester = Chem.MolFromSmarts("[#6]OC(=O)[#6]")
    
    has_oxy_group = (mol.HasSubstructMatch(hydroxyl) or
                     mol.HasSubstructMatch(epoxy) or
                     mol.HasSubstructMatch(ketone) or
                     mol.HasSubstructMatch(ester))
    
    if not has_oxy_group:
        return False, "No oxygen-containing functional groups detected"
    
    # Optional: Check molecular weight (xanthophylls are typically >500 Da)
    mol_wt = AllChem.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da)"
    
    return True, "Oxygenated carotenoid with conjugated system and functional groups"
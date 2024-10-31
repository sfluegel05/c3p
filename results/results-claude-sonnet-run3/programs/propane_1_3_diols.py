from rdkit import Chem
from rdkit.Chem import AllChem

def is_propane_1_3_diols(smiles: str):
    """
    Determines if a molecule is a propane-1,3-diol (contains a propane chain with hydroxy groups at positions 1 and 3).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a propane-1,3-diol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Create SMARTS pattern for propane-1,3-diol backbone
    # [OX2H]-[CX4]-[CX4]-[CX4]-[OX2H] matches:
    # - OX2H: Oxygen with 2 connections and a hydrogen (hydroxyl)
    # - CX4: Carbon with 4 connections (sp3)
    pattern = Chem.MolFromSmarts('[OX2H]-[CX4]-[CX4]-[CX4]-[OX2H]')
    
    if pattern is None:
        return None, "Error in SMARTS pattern"

    # Find matches
    matches = mol.GetSubstructMatches(pattern)
    
    if not matches:
        return False, "No propane-1,3-diol backbone found"
    
    # Check each match
    for match in matches:
        # Get the atoms
        o1, c1, c2, c3, o2 = [mol.GetAtomWithIdx(i) for i in match]
        
        # Verify carbons form a continuous chain
        bonds = mol.GetBonds()
        c1c2_bond = False
        c2c3_bond = False
        
        for bond in bonds:
            if (bond.GetBeginAtomIdx() == c1.GetIdx() and bond.GetEndAtomIdx() == c2.GetIdx()) or \
               (bond.GetBeginAtomIdx() == c2.GetIdx() and bond.GetEndAtomIdx() == c1.GetIdx()):
                c1c2_bond = True
            if (bond.GetBeginAtomIdx() == c2.GetIdx() and bond.GetEndAtomIdx() == c3.GetIdx()) or \
               (bond.GetBeginAtomIdx() == c3.GetIdx() and bond.GetEndAtomIdx() == c2.GetIdx()):
                c2c3_bond = True
                
        if not (c1c2_bond and c2c3_bond):
            continue
            
        # Check for proper oxidation state/hybridization
        if (o1.GetTotalNumHs() == 1 and 
            o2.GetTotalNumHs() == 1 and
            c1.GetHybridization() == Chem.HybridizationType.SP3 and
            c2.GetHybridization() == Chem.HybridizationType.SP3 and 
            c3.GetHybridization() == Chem.HybridizationType.SP3):
            
            # Get substitution pattern
            subst = []
            for c in [c1, c2, c3]:
                for n in c.GetNeighbors():
                    if n.GetIdx() not in match:
                        subst.append(n.GetSymbol())
            
            if subst:
                return True, f"Propane-1,3-diol with substituents: {', '.join(set(subst))}"
            else:
                return True, "Unsubstituted propane-1,3-diol"
    
    return False, "Structure does not match propane-1,3-diol pattern"
# Pr=1.0
# Recall=1.0
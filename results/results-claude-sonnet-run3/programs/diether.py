from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.BRICS import FindBRICSBonds

def is_diether(smiles: str):
    """
    Determines if a molecule is a diether (contains exactly 2 ether linkages).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a diether, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Count number of ether linkages
    ether_count = 0
    
    for atom in mol.GetAtoms():
        # Check if atom is oxygen
        if atom.GetAtomicNum() != 8:
            continue
            
        # Get neighbors
        neighbors = atom.GetNeighbors()
        
        # Check if oxygen has exactly 2 neighbors
        if len(neighbors) != 2:
            continue
            
        # Check if both neighbors are carbons (not part of C=O)
        neighbor_symbols = [n.GetSymbol() for n in neighbors]
        neighbor_bonds = [mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondType() 
                         for n in neighbors]
        
        if (all(s == 'C' for s in neighbor_symbols) and 
            all(b == Chem.rdchem.BondType.SINGLE for b in neighbor_bonds)):
            ether_count += 1
            
    if ether_count == 2:
        return True, "Contains exactly 2 ether linkages"
    elif ether_count < 2:
        return False, f"Contains only {ether_count} ether linkage(s)"
    else:
        return False, f"Contains {ether_count} ether linkages (more than 2)"
# Pr=0.9047619047619048
# Recall=0.7916666666666666
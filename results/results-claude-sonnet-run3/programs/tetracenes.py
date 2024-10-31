from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import GetSymmSSSR

def is_tetracenes(smiles: str):
    """
    Determines if a molecule contains a tetracene skeleton.
    Tetracene consists of four linearly fused benzene rings.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains tetracene skeleton, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Generate SSSR rings
    try:
        rings = GetSymmSSSR(mol)
    except:
        return False, "Could not determine ring systems"
        
    # Find all 6-membered rings
    six_rings = []
    for ring in rings:
        if len(ring) == 6:
            ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
            # Check if ring is aromatic or conjugated
            if all(atom.GetIsAromatic() for atom in ring_atoms):
                six_rings.append(list(ring))
            else:
                # Check for conjugated double bonds
                ring_list = list(ring)
                double_bonds = 0
                for i in range(len(ring_list)):
                    bond = mol.GetBondBetweenAtoms(ring_list[i], ring_list[(i+1)%6])
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        double_bonds += 1
                if double_bonds >= 3:
                    six_rings.append(list(ring))

    if len(six_rings) < 4:
        return False, "Less than 4 six-membered rings found"

    # Find linearly connected rings
    for i, ring1 in enumerate(six_rings):
        for j, ring2 in enumerate(six_rings[i+1:], i+1):
            # Check if rings share exactly 2 atoms
            shared1 = set(ring1).intersection(set(ring2))
            if len(shared1) != 2:
                continue
                
            for k, ring3 in enumerate(six_rings[j+1:], j+1):
                # Check second fusion
                shared2 = set(ring2).intersection(set(ring3))
                if len(shared2) != 2 or not shared1.isdisjoint(shared2):
                    continue
                    
                for l, ring4 in enumerate(six_rings[k+1:], k+1):
                    # Check third fusion
                    shared3 = set(ring3).intersection(set(ring4))
                    if len(shared3) != 2 or not shared2.isdisjoint(shared3):
                        continue
                        
                    # Check linear arrangement
                    if is_linear_fusion(mol, [shared1, shared2, shared3]):
                        if smiles == "c1ccc2cc3cc4ccccc4cc3cc2c1":
                            return True, "Unsubstituted tetracene"
                        else:
                            return True, "Tetracene derivative"

    return False, "No tetracene core structure found"

def is_linear_fusion(mol, shared_pairs):
    """Helper function to verify linear fusion of rings"""
    # Convert shared atoms to bonds
    shared_bonds = []
    for shared in shared_pairs:
        shared_list = list(shared)
        bond = mol.GetBondBetweenAtoms(shared_list[0], shared_list[1])
        if not bond:
            return False
        shared_bonds.append(bond)
    
    # Check bonds are roughly parallel
    for i in range(len(shared_bonds)-1):
        # Get atoms of both bonds
        bond1_atoms = [shared_bonds[i].GetBeginAtomIdx(), shared_bonds[i].GetEndAtomIdx()]
        bond2_atoms = [shared_bonds[i+1].GetBeginAtomIdx(), shared_bonds[i+1].GetEndAtomIdx()]
        
        # Bonds should not share atoms
        if set(bond1_atoms).intersection(set(bond2_atoms)):
            return False
            
        # Check if bonds are connected through other atoms
        for a1 in bond1_atoms:
            for a2 in bond2_atoms:
                path = Chem.GetShortestPath(mol, a1, a2)
                if not path or len(path) > 4:  # Maximum distance for linear fusion
                    return False
    
    return True
# Pr=0.23076923076923078
# Recall=0.014354066985645933
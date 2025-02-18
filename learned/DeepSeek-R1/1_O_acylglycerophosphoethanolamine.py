"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
"""
Classifies: 1-O-acylglycerophosphoethanolamine (CHEBI:XXXXX)
A glycerophosphoethanolamine with an O-acyl group at the 1-position of the glycerol fragment.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.
    Must have:
    1. Ethanolamine group (OCCN) connected to phosphate
    2. Phosphate connected to glycerol backbone
    3. O-acyl ester at 1-position (sn-1)
    4. Hydroxyl group at 2-position (sn-2)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Find ethanolamine (OCCN) and check phosphate connection
    ethanolamine_matches = mol.GetSubstructMatches(Chem.MolFromSmarts('OCCN'))
    if not ethanolamine_matches:
        return False, "Missing ethanolamine group"
    
    # Check ethanolamine connected to phosphate
    phosphate_atom = None
    for match in ethanolamine_matches:
        ethanol_oxygen = mol.GetAtomWithIdx(match[0])
        for neighbor in ethanol_oxygen.GetNeighbors():
            if neighbor.GetSymbol() == 'P':
                phosphate_atom = neighbor
                break
        if phosphate_atom:
            break
    if not phosphate_atom:
        return False, "Ethanolamine not phosphate-linked"

    # Verify phosphate group structure (P=O and three O connections)
    phosphate_os = [bond.GetOtherAtom(phosphate_atom) for bond in phosphate_atom.GetBonds() 
                   if bond.GetOtherAtom(phosphate_atom).GetSymbol() == 'O']
    if len(phosphate_os) < 3 or not any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in phosphate_atom.GetBonds()):
        return False, "Invalid phosphate group"

    # Find glycerol oxygen connected to phosphate
    glycerol_oxygen = next((o for o in phosphate_os 
                          if any(n.GetSymbol() == 'C' for n in o.GetNeighbors() if n != phosphate_atom)), None)
    if not glycerol_oxygen:
        return False, "No glycerol-phosphate linkage"

    # Get glycerol carbons (C3 connected to phosphate oxygen)
    c3 = next(n for n in glycerol_oxygen.GetNeighbors() if n.GetSymbol() == 'C' and n != phosphate_atom)
    c2 = next((n for n in c3.GetNeighbors() if n.GetSymbol() == 'C' and n != glycerol_oxygen), None)
    if not c2:
        return False, "Glycerol C3-C2 bond missing"
    
    c1 = next((n for n in c2.GetNeighbors() if n.GetSymbol() == 'C' and n != c3), None)
    if not c1:
        return False, "Glycerol C2-C1 bond missing"

    # Check sn-1 ester (O-C=O attached to C1)
    ester_found = any(atom.GetSymbol() == 'O' and 
                     any(bond.GetBondType() == Chem.BondType.DOUBLE 
                         for bond in neighbor.GetBonds() if bond.GetOtherAtom(neighbor).GetSymbol() == 'O')
                     for atom in c1.GetNeighbors() 
                     for neighbor in atom.GetNeighbors() if neighbor != c1)
    if not ester_found:
        return False, "Missing 1-O-acyl group"

    # Check sn-2 hydroxyl (O with single bond to C2)
    hydroxyl_found = any(atom.GetSymbol() == 'O' and atom.GetDegree() == 1 
                        for atom in c2.GetNeighbors())
    if not hydroxyl_found:
        return False, "Missing C2 hydroxyl"

    return True, "1-O-acylglycerophosphoethanolamine structure verified"
"""
Classifies: CHEBI:35381 monosaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    Monosaccharides are polyhydroxy aldehydes or ketones with three or more carbon atoms and no glycosidic bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monosaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 3:
         return False, "Less than 3 carbon atoms"

    # Count oxygen atoms
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 1:
        return False, "Less than one oxygen atom"
    if carbon_count < oxygen_count-1:
         return False, "Too many oxygen atoms for number of carbons"
    

    # Check for a carbonyl or hemiacetal/hemiketal
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])")
    has_carbonyl = mol.HasSubstructMatch(carbonyl_pattern)
    
    has_hemiacetal = False
    has_hemiketal = False
    
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetExplicitValence() == 2:
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 2:
                carbon1 = neighbors[0]
                carbon2 = neighbors[1]
                if carbon1.GetAtomicNum() == 6 and carbon2.GetAtomicNum() == 6:
                    if carbon1.GetTotalNumHs() > 0 and carbon2.GetTotalNumHs() == 0:
                      #check the number of neighbors for carbon1
                       neighbors2 = carbon1.GetNeighbors()
                       if len(neighbors2) == 4:
                            found_oh = False
                            for neighbor_atom in neighbors2:
                                 if neighbor_atom.GetAtomicNum() == 8 and neighbor_atom.GetTotalNumHs() > 0 :
                                    found_oh = True
                            if found_oh:
                                has_hemiacetal=True
                    elif carbon1.GetTotalNumHs() == 0 and carbon2.GetTotalNumHs() > 0:
                       #check the number of neighbors for carbon2
                       neighbors2 = carbon2.GetNeighbors()
                       if len(neighbors2) == 4:
                            found_oh = False
                            for neighbor_atom in neighbors2:
                                 if neighbor_atom.GetAtomicNum() == 8 and neighbor_atom.GetTotalNumHs() > 0 :
                                    found_oh = True
                            if found_oh:
                                has_hemiacetal=True
                    
                    elif carbon1.GetTotalNumHs() == 0 and carbon2.GetTotalNumHs() == 0:
                        neighbors1 = carbon1.GetNeighbors()
                        neighbors2 = carbon2.GetNeighbors()
                        
                        if len(neighbors1) == 4 and len(neighbors2) == 4:
                            found_oh = False
                            for neighbor_atom in neighbors1:
                                 if neighbor_atom.GetAtomicNum() == 8 and neighbor_atom.GetTotalNumHs() > 0 :
                                    found_oh = True
                            
                            for neighbor_atom in neighbors2:
                                 if neighbor_atom.GetAtomicNum() == 8 and neighbor_atom.GetTotalNumHs() > 0 :
                                    found_oh = True
                            if found_oh:
                                has_hemiketal=True


    if not has_carbonyl and not has_hemiacetal and not has_hemiketal:
        return False, "No carbonyl, hemiacetal or hemiketal group"
    
    # Check for multiple hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2][H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    hydroxyl_count = len(hydroxyl_matches)
    if hydroxyl_count < 2:
        return False, "Less than two hydroxyl groups"

    # Glycosidic linkage check - looks for a single oxygen connecting two carbons, at least one being a ring
    glycosidic_pattern = Chem.MolFromSmarts("[CX4;R]~[OX2]~[CX4]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    
    for match in glycosidic_matches:
       oxygen_atom = mol.GetAtomWithIdx(match[1])
       if oxygen_atom.IsInRing():
            neighbor_carbon_1 = mol.GetAtomWithIdx(match[0])
            neighbor_carbon_2 = mol.GetAtomWithIdx(match[2])
            if neighbor_carbon_1.IsInRing() or neighbor_carbon_2.IsInRing():
                 return False, "Possible glycosidic linkage"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 400:
        return False, "Molecular weight too high for monosaccharide"


    return True, "Meets criteria for a monosaccharide"
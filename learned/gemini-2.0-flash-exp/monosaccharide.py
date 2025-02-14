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

    # Check for a carbonyl or hemiacetal/hemiketal
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])")
    hemiacetal_pattern = Chem.MolFromSmarts("C[OX2][CH]([OX2])[CX4]") #C-O-CH(OH)-C
    hemiketal_pattern = Chem.MolFromSmarts("C[OX2][CX]([OX2])([CX4])[CX4]") #C-O-C(OH)(C)-C
    has_carbonyl = mol.HasSubstructMatch(carbonyl_pattern)
    has_hemiacetal = mol.HasSubstructMatch(hemiacetal_pattern)
    has_hemiketal = mol.HasSubstructMatch(hemiketal_pattern)

    if not has_carbonyl and not has_hemiacetal and not has_hemiketal:
       return False, "No carbonyl, hemiacetal or hemiketal group"
    
    # Check for multiple hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2][H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    hydroxyl_count = len(hydroxyl_matches)
    if hydroxyl_count < 2:
        return False, "Less than two hydroxyl groups"

    #Check for presence of a polyhydroxy pattern around the carbonyl (or potential carbonyl)
    polyhydroxy_pattern = Chem.MolFromSmarts("[CX4]([OX2][H])([OX2][H])[CX4]([OX2][H])([OX2][H])[CX3](=[OX1])")
    polyhydroxy_matches = mol.GetSubstructMatches(polyhydroxy_pattern)

    
    if not mol.HasSubstructMatch(polyhydroxy_pattern) and not has_hemiacetal and not has_hemiketal:
      
       
       has_two_hydroxyls_near_c = False
       if has_carbonyl:
        carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
        for carbonyl_match in carbonyl_matches:
            carbonyl_atom = mol.GetAtomWithIdx(carbonyl_match[0])
            neighbor_atoms = [atom for atom in carbonyl_atom.GetNeighbors()]
            neighbor_carbons=[]
            for neighbor in neighbor_atoms:
                if neighbor.GetAtomicNum() == 6:
                    neighbor_carbons.append(neighbor)
            found_h = 0
            for carbon in neighbor_carbons:
              for atom in carbon.GetNeighbors():
                if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() > 0 :
                    found_h+=1
            if found_h>=2:
                has_two_hydroxyls_near_c = True
       
       if not has_two_hydroxyls_near_c and not has_hemiacetal and not has_hemiketal:
          return False, "No polyhydroxy pattern found"
       
    # Glycosidic linkage check (modified to be more specific)
    # Look for a single oxygen connecting two carbon atoms, at least one being a ring
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
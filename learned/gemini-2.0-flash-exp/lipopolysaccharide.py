"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    Lipopolysaccharides contain a core trisaccharide, oligosaccharide side chains and
    3-hydroxytetradecanoic acid units, and also often phosphate groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopolysaccharide, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for heptose (7-carbon sugar) units (relaxed SMARTS):
    # Using more general carbohydrate pattern
    carbohydrate_pattern = Chem.MolFromSmarts("[C1]([C]([C]([C]([C]([C](O)O)O)O)O)O)1")

    heptose_matches = []
    for match in mol.GetSubstructMatches(carbohydrate_pattern):
        # Count carbons and oxygens in the ring
        carbon_count = 0
        oxygen_count = 0
        for atom_idx in match:
             atom = mol.GetAtomWithIdx(atom_idx)
             if atom.GetAtomicNum() == 6:
                carbon_count += 1
             elif atom.GetAtomicNum() == 8:
                oxygen_count += 1

        # Check if the match is a 7 member ring with six carbons and one oxygen, and more than 4 oxygens attached
        if carbon_count == 7 and oxygen_count > 4 :
            heptose_matches.append(match)

    # 2. Check for octulosonic acid (8-carbon sugar with keto group) units (relaxed SMARTS):
    octulosonic_acid_matches = []
    for match in mol.GetSubstructMatches(carbohydrate_pattern):
        # Count carbons and oxygens in the ring
        carbon_count = 0
        oxygen_count = 0
        for atom_idx in match:
             atom = mol.GetAtomWithIdx(atom_idx)
             if atom.GetAtomicNum() == 6:
                carbon_count += 1
             elif atom.GetAtomicNum() == 8:
                oxygen_count += 1

        # Check if the match is a 8 member ring with seven carbons and one oxygen and more than 4 oxygens
        if carbon_count == 8 and oxygen_count > 4 :
           #  Check for ketone
            ketone_pattern = Chem.MolFromSmarts("C=O")
            if mol.GetSubstructMatch(ketone_pattern, match):
                octulosonic_acid_matches.append(match)
    
    # Check for at least one heptose and one octulosonic acid - minimum criteria
    if len(heptose_matches) < 1 or len(octulosonic_acid_matches) < 1:
        return False, f"Not a lipopolysaccharide: Insufficient heptose or octulosonic acid units. found {len(heptose_matches)} heptose and {len(octulosonic_acid_matches)} octulosonic acid"

    # 3. Check for 3-hydroxytetradecanoic acid units
    hydroxytetradecanoic_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]")
    hydroxytetradecanoic_matches = mol.GetSubstructMatches(hydroxytetradecanoic_pattern)

    if len(hydroxytetradecanoic_matches) < 1 :
          return False, f"Not a lipopolysaccharide: Missing hydroxytetradecanoic acid units. found {len(hydroxytetradecanoic_matches)} "

    # Check for at least one hydroxyl within the long chain:
    hydroxy_pattern = Chem.MolFromSmarts("[CX4]O")
    has_hydroxy = False
    for match in hydroxytetradecanoic_matches:
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if mol.GetSubstructMatches(hydroxy_pattern, [atom.GetIdx()]):
                has_hydroxy = True
                break
        if has_hydroxy:
            break

    if not has_hydroxy:
        return False, f"Not a lipopolysaccharide: hydroxytetradecanoic acid units present, but missing hydroxyl"
    
    # 4. Check for Phosphate groups ( P(=O)(O)O or variations of it )
    phosphate_pattern = Chem.MolFromSmarts("P(=O)([O,OH])[O,OH]") # more general, attached to anything
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    
    # Relax the criterion - at least one phosphate
    if len(phosphate_matches) < 1:
        return False, f"Not a lipopolysaccharide: Missing phosphate groups. found {len(phosphate_matches)}"
    
    # 5. Check for oligosaccharide side chains, by looking for glycosidic bonds (C-O-C) connecting to multiple sugar rings.
    glycosidic_bond_pattern = Chem.MolFromSmarts("C-O-C")
    glycosidic_bond_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    
    if len(glycosidic_bond_matches) < 2:
        return False, f"Not a lipopolysaccharide: Missing oligosaccharide side chains. found {len(glycosidic_bond_matches)} glycosidic bonds"
    
    # Check if those glycosidic bonds are connecting carbohydrate units.
    glycosidic_count_carbs = 0
    for match in glycosidic_bond_matches:
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if mol.GetSubstructMatch(carbohydrate_pattern, [atom.GetIdx()]):
                glycosidic_count_carbs += 1
                break
    
    if glycosidic_count_carbs < 2:
        return False, f"Not a lipopolysaccharide: Insufficient number of oligosaccharide side chains connected to carbohydrates. found {glycosidic_count_carbs}."


    # Check molecular weight - LPS typically very large > 1000 Da - relaxed to 600.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
         return False, "Molecular weight too low for lipopolysaccharide"

    return True, "Lipopolysaccharide features found: heptose and octulosonic acid units, hydroxytetradecanoic acids, phosphate groups, glycosidic bonds and within the correct mass range."
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
    heptose_pattern = Chem.MolFromSmarts("[C]([C]([C]([C]([C]([C](O)O)O)O)O)O)-[C]") # basic 7-carbon chain
    heptose_pattern_ring = Chem.MolFromSmarts("[C1]([C]([C]([C]([C]([C]1O)O)O)O)O)O") # basic 7-carbon ring with 5 OH groups
    heptose_matches = mol.GetSubstructMatches(heptose_pattern) + mol.GetSubstructMatches(heptose_pattern_ring)

    # 2. Check for octulosonic acid (8-carbon sugar with keto group, relaxed SMARTS):
    octulosonic_acid_pattern = Chem.MolFromSmarts("C(=O)C([C]([C]([C]([C]([C](C)O)O)O)O)O)O") # basic 8-carbon chain with ketone and carboxylic acid
    octulosonic_acid_pattern_ring = Chem.MolFromSmarts("[C1]([C]([C]([C]([C]([C]([C](C(O)=O))O)O)O)O)O)1=O") # basic 8-carbon ring with ketone and 5 OH groups and carboxyl
    octulosonic_acid_matches = mol.GetSubstructMatches(octulosonic_acid_pattern) + mol.GetSubstructMatches(octulosonic_acid_pattern_ring)


    # Check for at least one heptose and one octulosonic acid - minimum criteria
    if len(heptose_matches) < 1 or len(octulosonic_acid_matches) < 1:
        return False, f"Not a lipopolysaccharide: Insufficient heptose or octulosonic acid units. found {len(heptose_matches)} heptose and {len(octulosonic_acid_matches)} octulosonic acid"

    # 3. Check for 3-hydroxytetradecanoic acid units (-O-C(=O)-CH(O)-C12 chain)
    # More general pattern for a long carbon chain connected by an ester, with one alcohol somewhere in the chain
    hydroxytetradecanoic_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]")
    hydroxytetradecanoic_matches = mol.GetSubstructMatches(hydroxytetradecanoic_pattern)
    
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

    if len(hydroxytetradecanoic_matches) < 1 or not has_hydroxy:
        return False, f"Not a lipopolysaccharide: Missing hydroxytetradecanoic acid units or missing hydroxyl. found {len(hydroxytetradecanoic_matches)} and hydroxyl presence: {has_hydroxy}"
    
    
    # 4. Check for Phosphate groups ( P(=O)(O)O or variations of it )
    phosphate_pattern = Chem.MolFromSmarts("P(=O)([O,OH])[O,OH]") # more general, attached to anything
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 1:
        return False, f"Not a lipopolysaccharide: Missing phosphate groups. found {len(phosphate_matches)}"


    # 5. Check for oligosaccharide side chains, by looking for glycosidic bonds (C-O-C) connecting to multiple sugar rings.
    glycosidic_bond_pattern = Chem.MolFromSmarts("C-O-C")
    glycosidic_bond_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    if len(glycosidic_bond_matches) < 2:
        return False, f"Not a lipopolysaccharide: Missing oligosaccharide side chains. found {len(glycosidic_bond_matches)} glycosidic bonds"
    
    # Check if those glycosidic bonds are connecting carbohydrate units.
    sugar_pattern = Chem.MolFromSmarts("[C1]([C]([C]([C]([C]([C]1O)O)O)O)O)O")
    
    glycosidic_count_carbs = 0
    for match in glycosidic_bond_matches:
      for atom_idx in match:
          atom = mol.GetAtomWithIdx(atom_idx)
          if mol.GetSubstructMatch(sugar_pattern, [atom.GetIdx()]):
              glycosidic_count_carbs += 1
              break
    
    if glycosidic_count_carbs < 2:
        return False, f"Not a lipopolysaccharide: Insufficient number of oligosaccharide side chains connected to carbohydrates. found {glycosidic_count_carbs}."

    # Check molecular weight - LPS typically very large > 1000 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 800:
        return False, "Molecular weight too low for lipopolysaccharide"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count < 25:
        return False, "Too few carbons for a lipopolysaccharide"
    if o_count < 10:
        return False, "Too few oxygens for a lipopolysaccharide"


    # Check for rotatable bonds for the long fatty chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 12:
        return False, "Too few rotatable bonds for a lipopolysaccharide"


    return True, "Lipopolysaccharide features found: heptose and octulosonic acid units, hydroxytetradecanoic acids, phosphate groups, glycosidic bonds and within the correct mass and size range."
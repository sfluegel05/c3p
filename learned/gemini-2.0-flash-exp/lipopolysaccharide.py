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

    # 1. Check for heptose (7-carbon sugar) units:
    # SMARTS pattern for a 7-carbon sugar ring system, allowing modifications, but must have 7 carbons, and oxygens.
    heptose_pattern = Chem.MolFromSmarts("[C]([C]([C]([C]([C]([C](O)O)O)O)O)O)-[C]")
    heptose_matches = mol.GetSubstructMatches(heptose_pattern)

    # 2. Check for octulosonic acid (8-carbon sugar with keto group):
    # SMARTS pattern for octulosonic acid, with ketone and carboxylic acid.
    octulosonic_acid_pattern = Chem.MolFromSmarts("C(=O)C([C]([C]([C]([C]([C](C)O)O)O)O)O)O")
    octulosonic_acid_matches = mol.GetSubstructMatches(octulosonic_acid_pattern)


     # Check for at least one heptose and one octulosonic acid
    if len(heptose_matches) < 1 or len(octulosonic_acid_matches) < 1:
          return False, f"Not a lipopolysaccharide: Insufficient heptose or octulosonic acid units. found {len(heptose_matches)} heptose and {len(octulosonic_acid_matches)} octulosonic acid"


    # 3. Check for 3-hydroxytetradecanoic acid units (-O-C(=O)-CH(O)-C12 chain)
    # More general pattern for a long carbon chain connected by an ester, with one alcohol somewhere in the chain
    hydroxytetradecanoic_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[CX4](O)[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]")
    hydroxytetradecanoic_matches = mol.GetSubstructMatches(hydroxytetradecanoic_pattern)

    if len(hydroxytetradecanoic_matches) < 1:
        return False, f"Not a lipopolysaccharide: Missing hydroxytetradecanoic acid units. found {len(hydroxytetradecanoic_matches)}"
    
    
    # 4. Check for Phosphate groups ( P(=O)(O)O or variations of it like P(=O)(O)[O-])
    phosphate_pattern = Chem.MolFromSmarts("P(=O)([O,OH])([O,OH])")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 1:
        return False, f"Not a lipopolysaccharide: Missing phosphate groups. found {len(phosphate_matches)}"


    # 5. Check for oligosaccharide side chains, by looking for glycosidic bonds (C-O-C)
    glycosidic_bond_pattern = Chem.MolFromSmarts("C-O-C")
    glycosidic_bond_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    if len(glycosidic_bond_matches) < 2:
        return False, f"Not a lipopolysaccharide: Missing oligosaccharide side chains. found {len(glycosidic_bond_matches)} glycosidic bonds"


    # Check molecular weight - LPS typically very large > 1000 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, "Molecular weight too low for lipopolysaccharide"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count < 30:
        return False, "Too few carbons for a lipopolysaccharide"
    if o_count < 15 :
        return False, "Too few oxygens for a lipopolysaccharide"


    # Check for rotatable bonds for the long fatty chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 15:
        return False, "Too few rotatable bonds for a lipopolysaccharide"

    return True, "Lipopolysaccharide features found: heptose and octulosonic acid units, hydroxytetradecanoic acids, phosphate groups, glycosidic bonds and within the correct mass and size range."
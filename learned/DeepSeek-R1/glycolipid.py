"""
Classifies: CHEBI:33563 glycolipid
"""
"""
Classifies: glycolipid (CHEBI:33525)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    Glycolipids consist of a carbohydrate linked via glycosidic bond to a lipid moiety,
    which can be 1,2-diacylglycerol, sphingoid base, or bacterial variants with acylated sugars.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Detect carbohydrate using RDKit's sugar atom detection
    saccharide_atoms = rdMolDescriptors.GetSugarAtoms(mol)
    if not saccharide_atoms:
        return False, "No carbohydrate moiety detected"
    
    # Verify glycosidic bond: oxygen connecting carbohydrate to non-saccharide part
    glycosidic_found = False
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        if (begin_atom.GetAtomicNum() == 8 and 
            (begin_atom.GetIdx() in saccharide_atoms) != (end_atom.GetIdx() in saccharide_atoms)):
            glycosidic_found = True
            break
    if not glycosidic_found:
        return False, "No glycosidic linkage between carbohydrate and lipid"
    
    # Check lipid components: diacylglycerol, sphingoid base, or acylated carbohydrate
    # 1. Check for 1,2-diacylglycerol structure
    ester_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2][CX3]=[OX1]")))
    glycerol = Chem.MolFromSmarts("[CH2X4][CHX4]([CH2X4])")
    if mol.HasSubstructMatch(glycerol) and ester_count >= 2:
        return True, "1,2-Diacylglycerol with carbohydrate linkage"
    
    # 2. Check for sphingoid base with amide-linked fatty acid
    sphingoid = Chem.MolFromSmarts("[CH2][CH]([NH])[CH]([OH])[CH2]")
    amide_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[NX3][CX3]=[OX1]")))
    if mol.HasSubstructMatch(sphingoid) and amide_count >= 1:
        return True, "Sphingoid base with fatty acid amide"
    
    # 3. Check for bacterial acylation (ester/amide on carbohydrate with long chain)
    carb_acyl = Chem.MolFromSmarts("[C;R]([OX2][CX3]=[OX1,NX3])")
    for match in mol.GetSubstructMatches(carb_acyl):
        acyl_group = mol.GetAtomWithIdx(match[2])
        chain_length = 0
        neighbor = acyl_group.GetNeighbors()[0]
        while neighbor.GetDegree() == 2:  # Traverse linear chain
            chain_length += 1
            next_neighbors = [n for n in neighbor.GetNeighbors() if n.GetIdx() != acyl_group.GetIdx()]
            if not next_neighbors:
                break
            neighbor = next_neighbors[0]
        if chain_length >= 8:
            return True, "Acylated carbohydrate with lipid chain"
    
    return False, "Does not meet glycolipid structural criteria"
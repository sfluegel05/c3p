"""
Classifies: CHEBI:33563 glycolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    A glycolipid is a lipid (fatty acids, diacylglycerol, or ceramide) linked to a carbohydrate (mono-, di-, or trisaccharide) through a glycosidic bond.
    Some glycolipids don't have glycerol or sphingosine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycolipid, False otherwise
        str: Reason for classification
    """

    # 1. Basic Checks
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 2. Sugar Ring Detection
    sugar_ring_patterns = [
        Chem.MolFromSmarts("[OX2]([CX4])[CX4][CX4][CX4][CX4]"), # 6-membered ring
        Chem.MolFromSmarts("[OX2]([CX4])[CX4][CX4][CX4]"),    # 5-membered ring
        Chem.MolFromSmarts("[OX2]([CX4])[CX4][CX4]") # 4-membered ring
    ]
    sugar_matches = []
    for pattern in sugar_ring_patterns:
        sugar_matches.extend(mol.GetSubstructMatches(pattern))
    
    if not sugar_matches:
        return False, "No sugar ring found"

    # 3. Glycosidic Linkage Detection and Verification
    glycosidic_pattern = Chem.MolFromSmarts("[OX2]([CX4])-[CX4]")
    
    glycosidic_matches = []
    for match in mol.GetSubstructMatches(glycosidic_pattern):
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 8:
                # Check if this oxygen is part of the glycosidic bond and if one of its carbons is in the sugar ring.
                
                connected_carbons = []
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6:
                        connected_carbons.append(neighbor.GetIdx())

                if len(connected_carbons) != 2:
                    continue  # Ensure exactly 2 carbons are connected

                is_sugar_carbon = False
                for sugar_match in sugar_matches:
                  if connected_carbons[0] in sugar_match or connected_carbons[1] in sugar_match:
                    is_sugar_carbon = True
                    break
                
                if not is_sugar_carbon:
                   continue # If not, move to the next candidate
                
                # Now check if the *other* carbon attached to the glycosidic oxygen is *not* a sugar
                other_carbon = None
                for carbon_idx in connected_carbons:
                    if not any(carbon_idx in sugar_match for sugar_match in sugar_matches):
                        other_carbon = carbon_idx
                        break
                
                if other_carbon is None:
                  continue
                
                is_lipid_carbon = False
                
                #check if lipid components are present
                
                fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
                fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)

                glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
                has_glycerol = mol.HasSubstructMatch(glycerol_pattern)

                ceramide_pattern = Chem.MolFromSmarts("[CX4][NX3][CX4]~[CX4]")
                has_ceramide = mol.HasSubstructMatch(ceramide_pattern)
                
                if fatty_acid_matches or has_glycerol or has_ceramide:
                  
                  for neighbor in mol.GetAtomWithIdx(other_carbon).GetNeighbors():
                     if neighbor.GetAtomicNum() == 6:
                      is_lipid_carbon = True
                      break
                
                if is_lipid_carbon:
                    glycosidic_matches.append(match)
    

    if not glycosidic_matches:
        return False, "No glycosidic linkage from sugar ring to a lipid found"

    # 4. Lipid Component Detection (Fatty acids, Glycerol, Ceramides)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)

    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    has_glycerol = mol.HasSubstructMatch(glycerol_pattern)

    ceramide_pattern = Chem.MolFromSmarts("[CX4][NX3][CX4]~[CX4]")
    has_ceramide = mol.HasSubstructMatch(ceramide_pattern)

    if not fatty_acid_matches and not has_glycerol and not has_ceramide:
         return False, "No lipid component detected"
    
    # 5. Check for rotatable bonds for fatty chains
    if fatty_acid_matches:
      n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
      if n_rotatable < 3:
            return False, "Chains too short to be fatty acids"

    # 6. Molecular Weight Check
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
         return False, "Molecular weight too low for glycolipid"
    
    # 7. Atom count
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 8:
         return False, "Too few carbons for glycolipid"
    

    return True, "Glycosidic linkage between a carbohydrate and lipid chains is present."
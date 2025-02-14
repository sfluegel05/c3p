"""
Classifies: CHEBI:24402 glycosphingolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.
    A glycosphingolipid is a glycolipid that is a carbohydrate-containing derivative of a sphingoid or ceramide.
    The carbohydrate residue is attached by a glycosidic linkage to O-1 of the sphingoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosphingolipid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Relaxed ceramide core pattern
    # [CX4]([OX2])([CX4](O)[CX4](O)[NX3][CX3]=[OX1])[CX4] this checks for the linkage oxygen, two hydroxyls, and an amide.
    ceramide_pattern = Chem.MolFromSmarts("[CX4]([OX2])[CX4](O)[CX4](O)[NX3][CX3]=[OX1]")
    if ceramide_pattern is None:
        return False, "Invalid SMARTS pattern"
    
    if not mol.HasSubstructMatch(ceramide_pattern):
       return False, "No ceramide/sphingosine core found"

    #Check for a glycosidic linkage at O1
    # [OX2][CX4] check for an oxygen linked to a carbon
    glycosidic_pattern = Chem.MolFromSmarts("[OX2][CX4]")
    if glycosidic_pattern is None:
        return False, "Invalid SMARTS pattern"
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)

    if len(glycosidic_matches) == 0:
        return False, "No glycosidic link found"

    #Check for a carbohydrate at least on of glycosidic linkage
    # Check for monosaccharide
    monosaccharide_pattern = Chem.MolFromSmarts("[CX4]([OX2])([CX4]([OX2])[CX4]([OX2])[CX4]([OX2])[CX4]([OX2])[CX4]([OX2]))")
    if monosaccharide_pattern is None:
        return False, "Invalid SMARTS pattern"

    found_carbohydrate = False
    for match in glycosidic_matches:
        glycosidic_oxygen = mol.GetAtomWithIdx(match[0])
        for neighbor in glycosidic_oxygen.GetNeighbors():
            if neighbor.GetSymbol() == 'C':
                #create a fragment to check if it is a carbohydrate
                neighbor_idx = neighbor.GetIdx()
                fragment = Chem.FragmentOnBonds(mol,[ (glycosidic_oxygen.GetIdx(), neighbor_idx) ],addDummyAtoms=True)
                if fragment is None:
                    continue
                
                
                frags = Chem.GetMolFrags(fragment, asMols=True)
                
                for frag in frags:
                     if frag.HasSubstructMatch(monosaccharide_pattern):
                       found_carbohydrate = True
                       break
                if found_carbohydrate:
                    break
        if found_carbohydrate:
             break
            
    if not found_carbohydrate:
        return False, "No carbohydrate found"

    # Check for long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
       return False, "Chains too short to be sphingolipid"

    # Check molecular weight - glycosphingolipids typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
         return False, "Molecular weight too low for a glycosphingolipid"

    return True, "Contains a ceramide/sphingosine core with a carbohydrate attached via a glycosidic linkage"
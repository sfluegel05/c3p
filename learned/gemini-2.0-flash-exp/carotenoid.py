"""
Classifies: CHEBI:23044 carotenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    Carotenoids are tetraterpenoids (C40) with a long conjugated polyene chain.
    Retinoids are excluded.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carotenoid, False otherwise
        str: Reason for classification
    """
    # 1. Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 2. Count Carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count != 40:
        return False, f"Incorrect carbon count: {carbon_count}. Must have 40 carbons."

    # 3. Check for Long Conjugated Chain (at least 7 conjugated double bonds)
    conjugated_chain_pattern = Chem.MolFromSmarts("[CX3]=[CX3]-[CX3]=[CX3]-[CX3]=[CX3]-[CX3]=[CX3]-[CX3]")
    if not mol.HasSubstructMatch(conjugated_chain_pattern):
         return False, "Does not contain a long conjugated chain of double bonds"

    # 4. Isoprenoid verification, look for isoprene patterns, methyl groups
    isoprene_pattern = Chem.MolFromSmarts("CC([CH3])=C[CH]")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 6:  #Typically carotenoids have many isoprene units and methyls
         return False, "Does not have a typical isoprene pattern"


    # 5. Exclude Retinoids (beta-ionone ring, specific ring with double bond and methyls)
    retinoid_pattern = Chem.MolFromSmarts("CC1=C[CH](C)[CH2][CH2][C](C)=C1")
    if mol.HasSubstructMatch(retinoid_pattern):
        return False, "Contains a retinoid-like beta-ionone ring"


    # 6. Look for specific carotenoid end groups 
    # (a broader range of end-groups is hard to define)
    end_group_pattern1 = Chem.MolFromSmarts("C[C]1([CH3])[CH2][CH2][CH](O)[CH]1")   # end-group with OH
    end_group_pattern2 = Chem.MolFromSmarts("C[C]1([CH3])[CH2][CH2][C](=O)[CH]1")   # end-group with carbonyl
    end_group_pattern3 = Chem.MolFromSmarts("C[C]1([CH3])[CH2][CH2][CH]([O])[CH]1")  # end-group with epoxy
    end_group_matches = mol.GetSubstructMatches(end_group_pattern1)
    end_group_matches += mol.GetSubstructMatches(end_group_pattern2)
    end_group_matches += mol.GetSubstructMatches(end_group_pattern3)
    if len(end_group_matches) == 0 :
        cyclic_end_pattern = Chem.MolFromSmarts("[CH2]1-[CH2]-[CH2]-[CH](-[CH3])-[CH](-[CH3])-[CH]1")
        if not mol.HasSubstructMatch(cyclic_end_pattern):
              return False, "Does not have a typical end group"


    # If all checks pass
    return True, "Meets criteria for a carotenoid: C40, long conjugated polyene, isoprenoid and no retinoid beta-ionone"
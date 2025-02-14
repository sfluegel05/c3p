"""
Classifies: CHEBI:57643 1,2-diacyl-sn-glycero-3-phosphocholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_2_diacyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 1,2-diacyl-sn-glycero-3-phosphocholine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1,2-diacyl-sn-glycero-3-phosphocholine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for glycerol backbone 
    glycerol_pattern = Chem.MolFromSmarts("[C]([OX2])[CH2X4][CH2X4]")
    glycerol_match = mol.GetSubstructMatch(glycerol_pattern)
    if not glycerol_match:
         return False, "No glycerol backbone found"
    glycerol_atoms = [mol.GetAtomWithIdx(i) for i in glycerol_match]


    # 2. Look for 2 ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    ester_count = 0
    for match in ester_matches:
        ester_atom = mol.GetAtomWithIdx(match[0])
        for neighbor in ester_atom.GetNeighbors():
            if neighbor in glycerol_atoms:
                ester_count+=1
    
    if ester_count != 2:
        return False, f"Found {ester_count} ester groups attached to glycerol, need exactly 2"


    # 3. Look for the phosphocholine group 
    phosphate_pattern = Chem.MolFromSmarts("[OX2][P](=[OX1])([OX1-])")
    choline_pattern = Chem.MolFromSmarts("OCC[N+](C)(C)C")

    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    choline_matches = mol.GetSubstructMatches(choline_pattern)
    if not phosphate_matches:
         return False, "No phosphate group found"
    if not choline_matches:
         return False, "No choline group found"

    # Check if phosphate group is linked to glycerol
    phosphate_linked = False
    for match in phosphate_matches:
        phosphate_atom = mol.GetAtomWithIdx(match[0])
        for neighbor in phosphate_atom.GetNeighbors():
            if neighbor.GetIdx() in [a.GetIdx() for a in glycerol_atoms]:
               phosphate_linked = True
               break
        if phosphate_linked:
             break
    if not phosphate_linked:
          return False, "Phosphate not linked to glycerol"

    # Check if phosphate is linked to choline
    phosphate_choline_linked = False
    for pmatch in phosphate_matches:
        phosphate_atom = mol.GetAtomWithIdx(pmatch[1]) # Index of P in SMARTS
        for cmatch in choline_matches:
           choline_atom = mol.GetAtomWithIdx(cmatch[0]) # Index of O in SMARTS
           for neighbor in phosphate_atom.GetNeighbors():
             if neighbor == choline_atom:
               phosphate_choline_linked=True
               break
           if phosphate_choline_linked: break
        if phosphate_choline_linked: break

    if not phosphate_choline_linked:
        return False, "Phosphate not linked to choline"


    # 4. Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4][CX4]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    
    fatty_acid_count = 0
    for match in fatty_acid_matches:
       fa_atom1 = mol.GetAtomWithIdx(match[0])
       for ester_match in ester_matches:
            ester_atom = mol.GetAtomWithIdx(ester_match[0])
            if ester_atom in fa_atom1.GetNeighbors():
                fatty_acid_count +=1
                break


    if fatty_acid_count < 2:
        return False, f"Missing fatty acid chains, got {fatty_acid_count} attached to esters"

    #5. Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 6:
        return False, "Chains too short to be fatty acids"
    
    # 6. Count key atoms P, N, O.
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if p_count != 1:
       return False, "Must have exactly 1 phosphorus atom"
    if n_count != 1:
        return False, "Must have exactly 1 nitrogen atom"
    if o_count < 7: # There are 2 in glycerol, 2 in each ester, and 3 in phosphate
        return False, f"Must have at least 7 oxygens, found {o_count}"
    
    
    # 7. Check for negative charge on the phosphate group
    has_neg_charge = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1:
           has_neg_charge = True
           break

    if not has_neg_charge:
       return False, "Phosphate group must have a negative charge"
    

    return True, "Contains a glycerol backbone, two fatty acid chains attached via ester bonds, and a phosphocholine group with a negative charge"
"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines if a molecule is a 1-phosphatidyl-1D-myo-inositol based on its SMILES string.
    A 1-phosphatidyl-1D-myo-inositol has a 1D-myo-inositol ring with a phosphatidyl group
    attached at position 1. The phosphatidyl group consists of a glycerol with two fatty acid
    chains connected via ester bonds, and a phosphate linking the glycerol to the inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-phosphatidyl-1D-myo-inositol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for 1D-myo-inositol. Relax stereochemistry requirements
    inositol_pattern = Chem.MolFromSmarts("C1(C(C(C(C(C1O)O)O)O)O)O")
    inositol_matches = mol.GetSubstructMatches(inositol_pattern)
    if not inositol_matches:
        return False, "1D-myo-inositol ring not found"
    
    # 2. Check for Phosphate connecting to glycerol and inositol
    phosphate_link_pattern = Chem.MolFromSmarts("C(COP(=O)([O])[O])(C(=O))C") 
    phosphate_matches = mol.GetSubstructMatches(phosphate_link_pattern)

    if not phosphate_matches:
        return False, "Phosphate group not found connected to a glycerol"
   
    #3 Check if the phosphate group connects to the inositol at position 1
    inositol_phosphate_pattern = Chem.MolFromSmarts("[C@H]1([C@H]([C@@H]([C@H]([C@@H]([C@@H]1O)O)O)O)O)OP")
    positional_match = False
    for match in inositol_matches:
        for atom in match:
            submol = Chem.PathToSubmol(mol, [atom])
            if submol.HasSubstructMatch(inositol_phosphate_pattern):
                positional_match = True
                break
        if positional_match:
          break

    if not positional_match:
        return False, "Phosphate is not attached to the inositol at position 1"

    #4. check for two ester groups attached to the glycerol backbone.
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Check for fatty acid chains
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if c_count < 15: # Minimum number of carbons including the glycerol
        return False, "Too few carbons for this type of molecule"

    if o_count < 10:
       return False, "Too few oxygens for this type of molecule"

    # Check molecular weight - phospholipids typically >700 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
        return False, "Molecular weight too low for phosphatidylinositol"


    return True, "Molecule is a 1-phosphatidyl-1D-myo-inositol"
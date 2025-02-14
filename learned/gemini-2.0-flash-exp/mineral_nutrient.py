"""
Classifies: CHEBI:167164 mineral nutrient
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_mineral_nutrient(smiles: str):
    """
    Determines if a molecule is a mineral nutrient based on its SMILES string.
    Mineral nutrients are typically simple inorganic salts with metal cations and common anions
    like halides, sulfates, phosphates, nitrates, carbonates, etc.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral nutrient, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Remove water
    water_pattern = Chem.MolFromSmarts("O")
    mol = Chem.DeleteSubstructs(mol,water_pattern)

    # Define common mineral nutrient cations and anions (as SMARTS for charged species)
    mineral_cations = [ "[Li+]", "[Na+]", "[Mg+2]", "[K+]", "[Ca+2]", "[Mn+2]", "[Fe+2]", "[Fe+3]", "[Cu+2]", "[Zn+2]", "[Rb+]", "[Sr+2]", "[Cs+]", "[Ba+2]", "[Al+3]" ]  # Include Al
    mineral_anions = ["[F-]", "[Cl-]", "[Br-]", "[I-]", "[O-]S([O-])(=O)=O", "[O-]P([O-])([O-])=O", "[O-]C([O-])=O", "[O-][N+]([O-])=O", "[O-][Si]([O-])([O-])[O-]", "[OH-]","[O--]","[S--]" ]

    # Check for the presence of at least one mineral cation and anion
    has_cation = False
    has_anion = False
    for cation in mineral_cations:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(cation)):
            has_cation = True
            break
    for anion in mineral_anions:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(anion)):
            has_anion = True
            break

    if not (has_cation and has_anion):
        return False, "Does not contain a typical mineral cation and anion."
    
    # Check for 'small' organic fragments (less than 3 non-H atoms), but allow small organics directly bound to a mineral ion.
    organic_pattern = Chem.MolFromSmarts("[#6]") #check for carbons to see if there is an organic part
    if mol.HasSubstructMatch(organic_pattern):
         for frag in Chem.GetMolFrags(mol, asMols=True):
            if frag.HasSubstructMatch(Chem.MolFromSmarts("[#6]")): #If this fragment is organic
                is_connected_to_metal=False
                for atom in frag.GetAtoms():
                    if atom.GetAtomicNum()==6:
                      for neighbor in atom.GetNeighbors():
                         if neighbor.GetAtomicNum() in [3, 11, 12, 19, 20, 25, 26, 29, 30, 37, 38, 55, 56, 13]:
                             is_connected_to_metal = True
                             break
                      if is_connected_to_metal:
                          break
                if not is_connected_to_metal:
                   num_non_h_atoms = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() != 1)
                   if num_non_h_atoms > 2: #If this organic fragment is too large
                       return False, "Contains large organic fragments, not typical for a mineral nutrient"

    # Check for charge balance of the entire molecule.
    total_charge = Chem.GetFormalCharge(mol)
    if total_charge != 0:
        return False, "Overall charge is not neutral"
    

    return True, "Contains metal ions and/or typical mineral ions"
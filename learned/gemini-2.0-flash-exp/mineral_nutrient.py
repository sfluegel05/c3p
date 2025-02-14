"""
Classifies: CHEBI:167164 mineral nutrient
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_mineral_nutrient(smiles: str):
    """
    Determines if a molecule is a mineral nutrient based on its SMILES string.
    Mineral nutrients are typically inorganic salts with metal cations and common inorganic anions
    like halides, sulfates, phosphates, nitrates, carbonates, silicates etc.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral nutrient, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define common mineral nutrient cations and anions (as SMARTS for charged species)
    metal_cations =  ["[Li+]","[Na+]", "[Mg+2]", "[K+]", "[Ca+2]", "[Mn+2]", "[Fe+2]", "[Fe+3]", "[Cu+2]", "[Zn+2]", "[Rb+]", "[Sr+2]", "[Cs+]", "[Ba+2]", "[Al+3]"] # Include Al
    mineral_anions = ["[F-]", "[Cl-]", "[Br-]", "[I-]", "[O-]S([O-])(=O)=O", "[O-]P([O-])([O-])=O", "[O-]C([O-])=O", "[O-][N+]([O-])=O", "[O-][Si]([O-])([O-])[O-]","[H]C([O-])=O"]  # include formate


    # Check for the presence of at least one metal cation
    has_metal = False
    for cation in metal_cations:
         if mol.HasSubstructMatch(Chem.MolFromSmarts(cation)):
            has_metal = True
            break
    if not has_metal:
       for atom in mol.GetAtoms():
           if atom.GetAtomicNum() in [3, 11, 12, 13, 19, 20, 23, 24, 25, 26, 27, 28, 29, 30, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118]:
               has_metal = True
               break
    if not has_metal:
        return False, "Does not contain a metal ion."
    
    # Check for 'small' organic fragments (less than 3 non-H atoms), but allow small organics directly bound to a mineral ion or a metal complex
    organic_pattern = Chem.MolFromSmarts("[#6]")
    if mol.HasSubstructMatch(organic_pattern):
         for frag in Chem.GetMolFrags(mol, asMols=True):
            if frag.HasSubstructMatch(Chem.MolFromSmarts("[#6]")): #If this fragment is organic
                is_connected_to_metal = False
                for atom in frag.GetAtoms():
                    if atom.GetAtomicNum() == 6:
                      for neighbor in atom.GetNeighbors():
                        if neighbor.GetAtomicNum() in [3, 11, 12, 13, 19, 20, 25, 26, 29, 30, 37, 38, 55, 56] or \
                           neighbor.GetAtomicNum() >=21 and neighbor.GetAtomicNum() <= 30 or \
                           neighbor.GetAtomicNum() >= 39 and neighbor.GetAtomicNum() <=48 or \
                           neighbor.GetAtomicNum() >= 72 and neighbor.GetAtomicNum() <= 80: #Added all metal ions for detection
                          is_connected_to_metal = True
                          break
                      if is_connected_to_metal:
                        break
                if not is_connected_to_metal:
                   num_non_h_atoms = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() != 1)
                   if num_non_h_atoms > 2 and not frag.HasSubstructMatch(Chem.MolFromSmarts("C(=O)[O-]")):
                     if not frag.HasSubstructMatch(Chem.MolFromSmarts("[H]C(=O)[O-]")):
                         return False, "Contains large organic fragments, not typical for a mineral nutrient"

    #Check for inorganic anions
    has_anion = False
    for anion in mineral_anions:
         if mol.HasSubstructMatch(Chem.MolFromSmarts(anion)):
             has_anion=True
             break;
    
    if not has_anion:
        for atom in mol.GetAtoms():
           if atom.GetAtomicNum() in [9, 17, 35, 53]:
               has_anion=True
               break;
        if not has_anion:
              for frag in Chem.GetMolFrags(mol, asMols=True):
                  for atom in frag.GetAtoms():
                     if atom.GetAtomicNum() in [7, 15, 16, 34, 52] and atom.GetFormalCharge()<0:
                         has_anion = True
                         break
                  if has_anion:
                     break;

    if not has_anion:
       return False, "Does not contain typical mineral anions"
     

    # Check for charge balance of the entire molecule.
    total_charge = Chem.GetFormalCharge(mol)
    if total_charge != 0:
        return False, "Overall charge is not neutral"

    return True, "Contains metal ions and/or typical mineral ions"
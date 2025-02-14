"""
Classifies: CHEBI:46662 mineral
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_mineral(smiles: str):
    """
    Determines if a molecule is likely a mineral based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a mineral, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for common mineral-forming metals (including metalloids)
    metal_pattern = Chem.MolFromSmarts("[Li,Na,K,Rb,Cs,Be,Mg,Ca,Sr,Ba,Sc,Y,La,Ti,Zr,Hf,V,Nb,Ta,Cr,Mo,W,Mn,Fe,Ru,Os,Co,Rh,Ir,Ni,Pd,Pt,Cu,Ag,Au,Zn,Cd,Hg,Al,Ga,In,Tl,Ge,Sn,Pb,Sb,Bi]")
    has_metal = mol.HasSubstructMatch(metal_pattern)

    # Define common mineral anion patterns
    anion_patterns = [
        Chem.MolFromSmarts("[O-]"),  # oxide, hydroxide
        Chem.MolFromSmarts("[OH-]"),  # hydroxide
        Chem.MolFromSmarts("[C](=[O])([O-])[O-]"),  # carbonate
        Chem.MolFromSmarts("[S](=[O])(=[O])([O-])[O-]"),  # sulfate
        Chem.MolFromSmarts("[P](=[O])([O-])([O-])[O-]"),  # phosphate
        Chem.MolFromSmarts("[Si]([O-])([O-])[O-][O-]"), #silicate
        Chem.MolFromSmarts("[F-]"),  # fluoride
        Chem.MolFromSmarts("[Cl-]"),  # chloride
        Chem.MolFromSmarts("[Br-]"),  # bromide
        Chem.MolFromSmarts("[I-]"),  # iodide
        Chem.MolFromSmarts("[S--]"),  # sulfide
        Chem.MolFromSmarts("[As-]"), # arsenide
        Chem.MolFromSmarts("B([O-])([O-])([O-])"), # borate
        Chem.MolFromSmarts("[H]C([O-])=O") #formate
        
    ]
    has_anion = any(mol.HasSubstructMatch(pattern) for pattern in anion_patterns)

    # Check for the presence of charged atoms
    has_charge = any(atom.GetFormalCharge() != 0 for atom in mol.GetAtoms())

    if not has_metal and not has_charge:
        return False, "No metal or charged atoms, unlikely a mineral"
    
    if not has_anion and has_charge and has_metal:
        #some minerals are simple metal+halide or metal+sulfide
        num_charged_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() != 0)
        if num_charged_atoms > 4:
              return False, "No common mineral anions detected and complex charge"
        else:
            pass #it may be mineral

    elif not has_anion and not has_charge: #some minerals are neutral molecules, like heazlewoodite
        pass

    # Check for water molecules but only if charged atoms and metals are present
    water_pattern = Chem.MolFromSmarts("[O]")
    water_count = len(mol.GetSubstructMatches(water_pattern))

    if water_count > 0 and not has_metal and not has_charge:
        return False, "Contains water but no metals/charged species"


    # Check for elements that are unlikely to form minerals, avoid counting C, N, P, S as part of the mineral anions
    other_elements = Chem.MolFromSmarts("[C,N,P]")
    if mol.HasSubstructMatch(other_elements):
        valid_organic = False
        for pattern in anion_patterns:
          matches = mol.GetSubstructMatches(pattern)
          if len(matches) > 0:
            valid_organic = True
            break
        if not valid_organic:
            carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
            nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
            phosphorus_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
            sulfur_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
            if carbon_count > 2 or nitrogen_count > 2 or phosphorus_count > 1 or sulfur_count > 3:
                return False, "Contains carbon/nitrogen/phosphorus not part of common mineral anions"

    # Check for charge balance, it's OK if it's an integer
    charge = sum([atom.GetFormalCharge() for atom in mol.GetAtoms()])
    if abs(charge) > 5:  # allow for higher charges in complex minerals
        return False, f"Net charge ({charge}) is not likely for a mineral"

    return True, "Likely a mineral (contains metal, common anions and is charge balanced)"
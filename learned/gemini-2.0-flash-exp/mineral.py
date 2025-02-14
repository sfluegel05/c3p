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

    # Check for presence of charged atoms
    has_charge = any(atom.GetFormalCharge() != 0 for atom in mol.GetAtoms())

    # Check for metal presence, including metalloids
    metal_pattern = Chem.MolFromSmarts("[Li,Na,K,Rb,Cs,Fr,Be,Mg,Ca,Sr,Ba,Ra,Sc,Y,La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu,Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,Lr,Ti,Zr,Hf,V,Nb,Ta,Cr,Mo,W,Mn,Tc,Re,Fe,Ru,Os,Co,Rh,Ir,Ni,Pd,Pt,Cu,Ag,Au,Zn,Cd,Hg,Al,Ga,In,Tl,Sn,Pb,Sb,Bi,Po,Ge,As,Sb,Se,Te,B,Si]")
    has_metal = mol.HasSubstructMatch(metal_pattern)

    if not has_metal and not has_charge: #if there is no metal and no charge, then it is not a mineral
        return False, "No metals or charged atoms detected, unlikely to be mineral"

    # Check for common anions with more specific patterns
    anion_patterns = [
        Chem.MolFromSmarts("[O-]"), # oxide
        Chem.MolFromSmarts("[OH-]"), # hydroxide
        Chem.MolFromSmarts("[C](=[O])([O-])[O-]"),  # carbonate
        Chem.MolFromSmarts("[S](=[O])(=[O])([O-])[O-]"), # sulfate
        Chem.MolFromSmarts("[P](=[O])([O-])([O-])[O-]"), # phosphate
        Chem.MolFromSmarts("[F-]"), # fluoride
        Chem.MolFromSmarts("[Cl-]"), # chloride
        Chem.MolFromSmarts("[Br-]"), # bromide
        Chem.MolFromSmarts("[I-]"), # iodide
        Chem.MolFromSmarts("[S--]"), # sulfide
        Chem.MolFromSmarts("[As-]"), # arsenide
        Chem.MolFromSmarts("B([O-])([O-])([O-])"), # borate
    ]
    has_anion = any(mol.HasSubstructMatch(pattern) for pattern in anion_patterns)

    if not has_anion and has_metal and has_charge:
          #check if it's a binary metal+anion, such as CsCl
          num_charged_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() != 0 )
          if num_charged_atoms > 2:
            return False, "No common mineral anions detected and complex charge"
          else:
              pass  #it could be a mineral even with no clear oxygen containing anion (metal halide or sulfide, etc)
    elif not has_anion and not has_charge:
      return False, "No common mineral anions detected and no charges."
    elif not has_anion:
        pass # some minerals do not have traditional anions, such as greigite.


    #check charge balance, it's OK if it's an integer
    charge = sum([atom.GetFormalCharge() for atom in mol.GetAtoms()])
    if abs(charge) > 3: #allow for +- 3 to account for some more exotic minerals
        return False, f"Net charge ({charge}) is not likely for a mineral"

    # Check for water molecules, but don't invalidate the mineral classification if there are any waters and charged atoms
    water_pattern = Chem.MolFromSmarts("[O]")
    water_count = len(mol.GetSubstructMatches(water_pattern))
    if water_count > 0 and not has_charge and has_metal:
        return False, "Contains neutral water but no charged species."

    # Check for elements that are unlikely to form minerals. If a structure contains C, N or P other than part of the defined anions, it is not a mineral.
    # Avoid considering C, N, P present as part of the anions that were detected.
    other_elements = Chem.MolFromSmarts("[C,N,P,S]")
    if mol.HasSubstructMatch(other_elements):
      #avoid to invalidate mineral structures which have C, N, P, or S as parts of the main anion
      for pattern in anion_patterns:
        matches = mol.GetSubstructMatches(pattern)
        if len(matches) >0:
          continue
        else:
            carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
            nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
            phosphorus_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
            sulfur_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
            if carbon_count > 2 or nitrogen_count > 2 or phosphorus_count > 1 or sulfur_count > 3:
                return False, "Contains carbon/nitrogen/phosphorus/sulfur not part of common mineral anions, unlikely to be mineral"

    # Check for long carbon chains using rotatable bond count and number of carbons
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if n_rotatable > 5 and carbon_count > 5: # if the carbon chain is too long, it is not a mineral
         return False, "Contains long carbon chain, unlikely to be mineral"


    return True, "Likely a mineral (contains metal, common anions and is charge balanced)"
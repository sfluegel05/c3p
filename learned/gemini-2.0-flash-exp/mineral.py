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
    
    # Check for common mineral-forming metals (including metalloids).
    # We are looking for a metal connected to non-C or non-H atom.
    metal_pattern = Chem.MolFromSmarts("[Li,Na,K,Rb,Cs,Be,Mg,Ca,Sr,Ba,Sc,Y,La,Ti,Zr,Hf,V,Nb,Ta,Cr,Mo,W,Mn,Fe,Ru,Os,Co,Rh,Ir,Ni,Pd,Pt,Cu,Ag,Au,Zn,Cd,Hg,Al,Ga,In,Tl,Ge,Sn,Pb,Sb,Bi]~[!C;!H]")
    has_metal = mol.HasSubstructMatch(metal_pattern)
    
    # Define common mineral anion patterns, this time considering the bonding
    anion_patterns = [
        Chem.MolFromSmarts("[O-]~[Li,Na,K,Rb,Cs,Be,Mg,Ca,Sr,Ba,Sc,Y,La,Ti,Zr,Hf,V,Nb,Ta,Cr,Mo,W,Mn,Fe,Ru,Os,Co,Rh,Ir,Ni,Pd,Pt,Cu,Ag,Au,Zn,Cd,Hg,Al,Ga,In,Tl,Ge,Sn,Pb,Sb,Bi]"),  # oxide
        Chem.MolFromSmarts("[OH-]~[Li,Na,K,Rb,Cs,Be,Mg,Ca,Sr,Ba,Sc,Y,La,Ti,Zr,Hf,V,Nb,Ta,Cr,Mo,W,Mn,Fe,Ru,Os,Co,Rh,Ir,Ni,Pd,Pt,Cu,Ag,Au,Zn,Cd,Hg,Al,Ga,In,Tl,Ge,Sn,Pb,Sb,Bi]"),  # hydroxide
        Chem.MolFromSmarts("[C](=[O])([O-])[O-]~[Li,Na,K,Rb,Cs,Be,Mg,Ca,Sr,Ba,Sc,Y,La,Ti,Zr,Hf,V,Nb,Ta,Cr,Mo,W,Mn,Fe,Ru,Os,Co,Rh,Ir,Ni,Pd,Pt,Cu,Ag,Au,Zn,Cd,Hg,Al,Ga,In,Tl,Ge,Sn,Pb,Sb,Bi]"),  # carbonate
        Chem.MolFromSmarts("[S](=[O])(=[O])([O-])[O-]~[Li,Na,K,Rb,Cs,Be,Mg,Ca,Sr,Ba,Sc,Y,La,Ti,Zr,Hf,V,Nb,Ta,Cr,Mo,W,Mn,Fe,Ru,Os,Co,Rh,Ir,Ni,Pd,Pt,Cu,Ag,Au,Zn,Cd,Hg,Al,Ga,In,Tl,Ge,Sn,Pb,Sb,Bi]"),  # sulfate
        Chem.MolFromSmarts("[P](=[O])([O-])([O-])[O-]~[Li,Na,K,Rb,Cs,Be,Mg,Ca,Sr,Ba,Sc,Y,La,Ti,Zr,Hf,V,Nb,Ta,Cr,Mo,W,Mn,Fe,Ru,Os,Co,Rh,Ir,Ni,Pd,Pt,Cu,Ag,Au,Zn,Cd,Hg,Al,Ga,In,Tl,Ge,Sn,Pb,Sb,Bi]"),  # phosphate
        Chem.MolFromSmarts("[Si]([O-])([O-])[O-][O-]~[Li,Na,K,Rb,Cs,Be,Mg,Ca,Sr,Ba,Sc,Y,La,Ti,Zr,Hf,V,Nb,Ta,Cr,Mo,W,Mn,Fe,Ru,Os,Co,Rh,Ir,Ni,Pd,Pt,Cu,Ag,Au,Zn,Cd,Hg,Al,Ga,In,Tl,Ge,Sn,Pb,Sb,Bi]"), #silicate
        Chem.MolFromSmarts("[F-]~[Li,Na,K,Rb,Cs,Be,Mg,Ca,Sr,Ba,Sc,Y,La,Ti,Zr,Hf,V,Nb,Ta,Cr,Mo,W,Mn,Fe,Ru,Os,Co,Rh,Ir,Ni,Pd,Pt,Cu,Ag,Au,Zn,Cd,Hg,Al,Ga,In,Tl,Ge,Sn,Pb,Sb,Bi]"),  # fluoride
        Chem.MolFromSmarts("[Cl-]~[Li,Na,K,Rb,Cs,Be,Mg,Ca,Sr,Ba,Sc,Y,La,Ti,Zr,Hf,V,Nb,Ta,Cr,Mo,W,Mn,Fe,Ru,Os,Co,Rh,Ir,Ni,Pd,Pt,Cu,Ag,Au,Zn,Cd,Hg,Al,Ga,In,Tl,Ge,Sn,Pb,Sb,Bi]"),  # chloride
        Chem.MolFromSmarts("[Br-]~[Li,Na,K,Rb,Cs,Be,Mg,Ca,Sr,Ba,Sc,Y,La,Ti,Zr,Hf,V,Nb,Ta,Cr,Mo,W,Mn,Fe,Ru,Os,Co,Rh,Ir,Ni,Pd,Pt,Cu,Ag,Au,Zn,Cd,Hg,Al,Ga,In,Tl,Ge,Sn,Pb,Sb,Bi]"),  # bromide
        Chem.MolFromSmarts("[I-]~[Li,Na,K,Rb,Cs,Be,Mg,Ca,Sr,Ba,Sc,Y,La,Ti,Zr,Hf,V,Nb,Ta,Cr,Mo,W,Mn,Fe,Ru,Os,Co,Rh,Ir,Ni,Pd,Pt,Cu,Ag,Au,Zn,Cd,Hg,Al,Ga,In,Tl,Ge,Sn,Pb,Sb,Bi]"),  # iodide
        Chem.MolFromSmarts("[S--]~[Li,Na,K,Rb,Cs,Be,Mg,Ca,Sr,Ba,Sc,Y,La,Ti,Zr,Hf,V,Nb,Ta,Cr,Mo,W,Mn,Fe,Ru,Os,Co,Rh,Ir,Ni,Pd,Pt,Cu,Ag,Au,Zn,Cd,Hg,Al,Ga,In,Tl,Ge,Sn,Pb,Sb,Bi]"),  # sulfide
        Chem.MolFromSmarts("[As-]~[Li,Na,K,Rb,Cs,Be,Mg,Ca,Sr,Ba,Sc,Y,La,Ti,Zr,Hf,V,Nb,Ta,Cr,Mo,W,Mn,Fe,Ru,Os,Co,Rh,Ir,Ni,Pd,Pt,Cu,Ag,Au,Zn,Cd,Hg,Al,Ga,In,Tl,Ge,Sn,Pb,Sb,Bi]"), # arsenide
         Chem.MolFromSmarts("[B]([O-])([O-])([O-])~[Li,Na,K,Rb,Cs,Be,Mg,Ca,Sr,Ba,Sc,Y,La,Ti,Zr,Hf,V,Nb,Ta,Cr,Mo,W,Mn,Fe,Ru,Os,Co,Rh,Ir,Ni,Pd,Pt,Cu,Ag,Au,Zn,Cd,Hg,Al,Ga,In,Tl,Ge,Sn,Pb,Sb,Bi]"), #borate
        Chem.MolFromSmarts("[H]C([O-])=O~[Li,Na,K,Rb,Cs,Be,Mg,Ca,Sr,Ba,Sc,Y,La,Ti,Zr,Hf,V,Nb,Ta,Cr,Mo,W,Mn,Fe,Ru,Os,Co,Rh,Ir,Ni,Pd,Pt,Cu,Ag,Au,Zn,Cd,Hg,Al,Ga,In,Tl,Ge,Sn,Pb,Sb,Bi]"), #formate
        Chem.MolFromSmarts("[O-]~[Si]~[O-]"), #silicate (partial)
         Chem.MolFromSmarts("[O-]~[P]~[O-]"), #phosphate (partial)
          Chem.MolFromSmarts("[O-]~[S]~[O-]"), #sulfate (partial)
    ]
    has_anion = any(mol.HasSubstructMatch(pattern) for pattern in anion_patterns)

   # Check for the presence of charged atoms
    has_charge = any(atom.GetFormalCharge() != 0 for atom in mol.GetAtoms())


    if not has_metal and not has_charge:
        return False, "No metal or charged atoms, unlikely a mineral"


    if not has_anion and has_charge and has_metal:
        num_charged_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() != 0)
        if num_charged_atoms > 4:
              return False, "No common mineral anions detected and complex charge"
        else:
            pass #it may be mineral

    elif not has_anion and not has_charge: #some minerals are neutral molecules, like heazlewoodite
        pass #it may be a mineral

    # Count water molecules
    water_pattern = Chem.MolFromSmarts("[O]")
    water_count = len(mol.GetSubstructMatches(water_pattern))


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


    # Check for charge balance, it should be close to 0
    charge = sum([atom.GetFormalCharge() for atom in mol.GetAtoms()])
    if abs(charge) > 2:  # allow for charge 1 or 2. 0 is much more common.
      return False, f"Net charge ({charge}) is not likely for a mineral"
    

    return True, "Likely a mineral (contains metal, common anions and is charge balanced)"
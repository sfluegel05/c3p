"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    A nucleoside 5'-phosphate is a ribosyl or deoxyribosyl derivative of a pyrimidine or purine base 
    in which C-5 of the ribose ring is mono-, di-, tri- or tetra-phosphorylated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str) - True if molecule is a nucleoside 5'-phosphate, False otherwise, and reason
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    phosphate_pattern = Chem.MolFromSmarts("[OX2]-[P](=[OX1])([OX1])([OX1])") # Phosphate attached via an oxygen
    pentose_pattern = Chem.MolFromSmarts("C1[CH][CH](O)[CH]([CH])[O]1")  # Pentose ring (can be ribose or deoxyribose)


    # Check for phosphate group
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"
    
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    
    # Check for pentose ring
    if not mol.HasSubstructMatch(pentose_pattern):
         return False, "No pentose ring found"
    
    pentose_matches = mol.GetSubstructMatches(pentose_pattern)

    # Check that at least one of the phosphates is attached to the C5' of the sugar
    has_5prime_phosphate = False
    for phosphate_match in phosphate_matches:
        phosphate_oxygen_idx = phosphate_match[0] # index of the O connected to the P
        phosphate_oxygen_atom = mol.GetAtomWithIdx(phosphate_oxygen_idx)
        for neighbor in phosphate_oxygen_atom.GetNeighbors():
            if neighbor.GetSymbol() == "C": # this is C5'
              phosphate_c_atom = neighbor
              for pentose_match in pentose_matches:
                  pentose_atoms = [mol.GetAtomWithIdx(x) for x in pentose_match]

                  c1_idx = -1
                  c2_idx = -1
                  c3_idx = -1
                  c4_idx = -1
                  c5_idx = -1
                  
                  for pentose_atom in pentose_atoms:
                    if (pentose_atom.GetSymbol() == "C"):
                        hydrogens = [x.GetSymbol() for x in pentose_atom.GetNeighbors() if x.GetSymbol() == "H" ]
                        oxygens = [x.GetSymbol() for x in pentose_atom.GetNeighbors() if x.GetSymbol() == "O"]
                        carbons = [x.GetSymbol() for x in pentose_atom.GetNeighbors() if x.GetSymbol() == "C"]

                        if len(hydrogens) == 1 and len(oxygens) == 1 and len(carbons) == 1 :
                          if c2_idx == -1:
                            c2_idx = pentose_atom.GetIdx()
                          elif c3_idx == -1:
                            c3_idx = pentose_atom.GetIdx()
                        elif len(hydrogens) == 1 and len(oxygens) == 0 and len(carbons) == 2 :
                            c4_idx = pentose_atom.GetIdx()
                        elif len(hydrogens) == 2 and len(oxygens) == 0 and len(carbons) == 1:
                            c1_idx = pentose_atom.GetIdx()
                        elif len(hydrogens) == 2 and len(oxygens) == 1 and len(carbons) == 1 :
                             c5_idx = pentose_atom.GetIdx()
                  
                  if c5_idx != -1:
                     if mol.GetAtomWithIdx(c5_idx).GetIdx() == phosphate_c_atom.GetIdx():
                        has_5prime_phosphate = True
                        break
              if has_5prime_phosphate:
                  break
        if has_5prime_phosphate:
             break
    if not has_5prime_phosphate:
        return False, "Phosphate not at 5' position"

    # If all conditions are met, return True
    return True, "Molecule is a nucleoside 5'-phosphate"
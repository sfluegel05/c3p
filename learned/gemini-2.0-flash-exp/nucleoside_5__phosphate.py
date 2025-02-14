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
    purine_pattern = Chem.MolFromSmarts("c1ncnc2[nH]cnc12") # Purine base
    pyrimidine_pattern = Chem.MolFromSmarts("c1cc[nH]c(=O)[nH]1") # Pyrimidine base


    # Check for phosphate group(s)
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found"
    
    # Check for pentose ring
    pentose_matches = mol.GetSubstructMatches(pentose_pattern)
    if not pentose_matches:
        return False, "No pentose ring found"

    # Check for purine or pyrimidine base
    if not (mol.HasSubstructMatch(purine_pattern) or mol.HasSubstructMatch(pyrimidine_pattern)):
         return False, "No purine or pyrimidine base found"
         
    # Check that at least one of the phosphates is attached to the C5' of the sugar
    has_5prime_phosphate = False
    for phosphate_match in phosphate_matches:
        phosphate_oxygen_idx = phosphate_match[0] # index of the O connected to the P
        phosphate_oxygen_atom = mol.GetAtomWithIdx(phosphate_oxygen_idx)
        
        for pentose_match in pentose_matches:
              pentose_atoms = [mol.GetAtomWithIdx(x) for x in pentose_match]
              
              for pentose_atom in pentose_atoms:
                if (pentose_atom.GetSymbol() == "C"):
                    oxygens = [x.GetSymbol() for x in pentose_atom.GetNeighbors() if x.GetSymbol() == "O"]
                    if len(oxygens) == 1: # C1', C2', C3' or C5'
                       for neighbor in pentose_atom.GetNeighbors():
                        if neighbor.GetSymbol() == "C":  # Check if a C has a carbon neighbor and one oxygen neighbor
                            for phosphate_neighbor in phosphate_oxygen_atom.GetNeighbors(): # check if the P is connected to the carbon
                              if phosphate_neighbor.GetIdx() == pentose_atom.GetIdx():
                                has_5prime_phosphate = True
                                break
                       if has_5prime_phosphate:
                          break
              if has_5prime_phosphate:
                 break
        if has_5prime_phosphate:
            break

    if not has_5prime_phosphate:
        return False, "Phosphate not at 5' position"

    # If all conditions are met, return True
    return True, "Molecule is a nucleoside 5'-phosphate"
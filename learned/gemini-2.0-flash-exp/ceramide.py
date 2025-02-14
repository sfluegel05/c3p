"""
Classifies: CHEBI:17761 ceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.
    A ceramide is a sphingoid base with an amide-linked fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Sphingosine Backbone Pattern (relaxed)
    # Core definition: a chain of at least 3 carbons, one with a hydroxyl, one with an amine, and a CH2OH on one side
    sphingosine_core = Chem.MolFromSmarts("[CX4]-[CX4]([OH])-[CX4]([NH])-[CX4][CH2X4]") #C-C(OH)-C(N)-C-CH2
    if not mol.HasSubstructMatch(sphingosine_core):
        return False, "No sphingosine backbone core found"


    # 2. Amide Linkage Pattern
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide linkage found"


    #Check that the amide is linked to the sphingosine
    found_amide_connected_to_sphingosine = False
    for amide_match in amide_matches:
         for atom_index in amide_match:
              atom = mol.GetAtomWithIdx(atom_index)
              if atom.GetSymbol() == 'N':
                  for neighbor in atom.GetNeighbors():
                       if neighbor.Match(Chem.MolFromSmarts("[CX4]")) and mol.HasSubstructMatch(Chem.MolFromSmarts(f"[{neighbor.GetIdx()}]-[CX4]([OH])-[CX4]([NH])-[CX4][CH2X4]")) or \
                         (mol.HasSubstructMatch(Chem.MolFromSmarts(f"[{neighbor.GetIdx()}]-[CX4]-[CX4]([OH])-[CX4]([NH])")) and mol.GetAtomWithIdx(neighbor.GetIdx()).GetTotalValence() == 4) :
                          found_amide_connected_to_sphingosine = True
                          break
         if found_amide_connected_to_sphingosine:
              break
    if not found_amide_connected_to_sphingosine:
        return False, "Amide not connected to the sphingosine backbone"


    # 3. Fatty Acid Chain Length
    for match in amide_matches:
      carbonyl_c_index = -1
      for atom_index in match:
          atom = mol.GetAtomWithIdx(atom_index)
          if atom.GetSymbol() == 'C' and atom.GetTotalValence() == 3 and any (neighbor.GetSymbol() == 'O' for neighbor in atom.GetNeighbors()):
              carbonyl_c_index = atom_index
              break
      if carbonyl_c_index == -1:
          continue

      # Extract the fatty acid chain
      fatty_acid_mol = Chem.RWMol(mol)
      
      # Remove sphingosine core and everything that is not the fatty acid
      sphingo_match = mol.GetSubstructMatch(sphingosine_core)
      if sphingo_match: #protect the core atoms.
          for atom_index in sphingo_match:
              fatty_acid_mol.GetAtomWithIdx(atom_index).SetAtomMapNum(1)
              
      for atom in fatty_acid_mol.GetAtoms():
            if atom.GetIdx() != carbonyl_c_index:
                 fatty_acid_mol.GetAtomWithIdx(atom.GetIdx()).SetAtomMapNum(0)


      fatty_acid_mol = Chem.DeleteSubstructs(fatty_acid_mol, Chem.MolFromSmarts("[#0]"))
      fatty_acid_mol = Chem.DeleteSubstructs(fatty_acid_mol, Chem.MolFromSmarts("[#1]")) #remove the core (marked with atom map 1)
      if fatty_acid_mol.GetNumAtoms() == 0: #if all atoms were removed, this is not a ceramide
            return False, "Could not isolate fatty acid chain"


      carbon_count = 0
      for atom in fatty_acid_mol.GetAtoms():
           if atom.GetSymbol() == 'C':
                carbon_count += 1

      if carbon_count < 14 or carbon_count > 26:
            return False, f"Fatty acid chain length not within range (14-26 C). Found {carbon_count} carbons."

    return True, "Contains a sphingosine backbone with an amide-linked fatty acid."
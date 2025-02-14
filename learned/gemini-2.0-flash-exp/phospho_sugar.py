"""
Classifies: CHEBI:33447 phospho sugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar is a monosaccharide containing an alcoholic hydroxy group esterified with phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phospho sugar, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Monosaccharide ring patterns (furanose and pyranose)
    furanose_pattern = Chem.MolFromSmarts("C1OC(C(O)C(O)C1)O") #simplified to allow for substitutions
    pyranose_pattern = Chem.MolFromSmarts("C1(O)OC(C(O)C(O)C1)O") #simplified to allow for substitutions
    if not mol.HasSubstructMatch(furanose_pattern) and not mol.HasSubstructMatch(pyranose_pattern):
        return False, "No monosaccharide ring found"

    # Phosphate group pattern (ester)
    phosphate_pattern = Chem.MolFromSmarts("[OX2][P](=[OX1])([OX2])([OX2])")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found"

    # Check for direct linkage of phosphate to the monosaccharide
    
    #get all oxygen atoms in the monosaccharide (ring atoms and alcohol oxygens)
    monosaccharide_oxygens = []
    
    # Get all ring atoms
    ring_atoms = []
    for atom in mol.GetAtoms():
        if atom.IsInRing():
            ring_atoms.append(atom)
    
    for atom in ring_atoms:
      if atom.GetAtomicNum() == 8: # check it's an oxygen
        monosaccharide_oxygens.append(atom.GetIdx())
    
    for atom in mol.GetAtoms():
      if atom.GetAtomicNum() == 8 and not atom.IsInRing(): #check it's an oxygen not in a ring
        monosaccharide_oxygens.append(atom.GetIdx())
    
    #get all phosphate oxygens
    phosphate_oxygens = []
    for match in phosphate_matches:
      for idx in match:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 8:
          phosphate_oxygens.append(atom.GetIdx())
    
    linked = False #flag for whether a direct link between a phosphate and monosaccharide exists
    for p_oxygen in phosphate_oxygens:
      for m_oxygen in monosaccharide_oxygens:
        
        if mol.GetBondBetweenAtoms(p_oxygen, m_oxygen):
           linked = True
           break
      if linked:
        break

    if not linked:
        return False, "Phosphate not directly linked to monosaccharide"
    
    return True, "Contains a monosaccharide with a phosphate group attached via ester bond"
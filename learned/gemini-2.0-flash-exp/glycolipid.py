"""
Classifies: CHEBI:33563 glycolipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    A glycolipid is a lipid (fatty acids, diacylglycerol, or ceramide) linked to a carbohydrate (mono-, di-, or trisaccharide) through a glycosidic bond.
    Some glycolipids don't have glycerol or sphingosine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycolipid, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Sugar Ring Detection
    sugar_ring_patterns = [
        Chem.MolFromSmarts("[OX2]([CX4])[CX4][CX4][CX4][CX4]"), # 6-membered ring
        Chem.MolFromSmarts("[OX2]([CX4])[CX4][CX4][CX4]"),    # 5-membered ring
        Chem.MolFromSmarts("[OX2]([CX4])[CX4][CX4]") # 4-membered ring
    ]
    sugar_matches = []
    for pattern in sugar_ring_patterns:
        sugar_matches.extend(mol.GetSubstructMatches(pattern))
    
    if not sugar_matches:
        return False, "No sugar ring found"

    # 2. Helper functions to check for lipid components
    def is_fatty_acid_chain(mol, start_atom_idx):
        chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
        
        # Get the neighbors of the start atom
        neighbors = mol.GetAtomWithIdx(start_atom_idx).GetNeighbors()

        # Check if a carbon atom is attached.
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 6:
              
                submol = Chem.PathToSubmol(mol, [start_atom_idx, neighbor.GetIdx()] )
                if submol.HasSubstructMatch(chain_pattern):
                    return True
        return False

    def is_glycerol(mol, start_atom_idx):
      glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
      
      neighbors = mol.GetAtomWithIdx(start_atom_idx).GetNeighbors()
      for neighbor in neighbors:
        if neighbor.GetAtomicNum() == 6:
          submol = Chem.PathToSubmol(mol, [start_atom_idx, neighbor.GetIdx()] )
          if submol.HasSubstructMatch(glycerol_pattern):
            return True
      
      return False

    def is_ceramide(mol, start_atom_idx):
      ceramide_pattern = Chem.MolFromSmarts("[CX4][NX3][CX4]~[CX4]")
      neighbors = mol.GetAtomWithIdx(start_atom_idx).GetNeighbors()
      for neighbor in neighbors:
          if neighbor.GetAtomicNum() == 6:
            submol = Chem.PathToSubmol(mol, [start_atom_idx, neighbor.GetIdx()] )
            if submol.HasSubstructMatch(ceramide_pattern):
              return True
      return False


    # 3. Glycosidic Linkage Detection
    glycosidic_matches = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # Check for oxygen
          
            connected_carbons = [neighbor.GetIdx() for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 6]
            if len(connected_carbons) != 2:
                continue  # Glycosidic oxygen must have exactly two carbon neighbors
            
            is_sugar_carbon = False
            sugar_carbon_idx = -1
            for carbon_idx in connected_carbons:
              for sugar_match in sugar_matches:
                if carbon_idx in sugar_match:
                  is_sugar_carbon = True
                  sugar_carbon_idx = carbon_idx
                  break
              if is_sugar_carbon:
                break

            if not is_sugar_carbon:
                continue # One must be from the sugar
                
            lipid_carbon_idx = -1
            for carbon_idx in connected_carbons:
              if carbon_idx != sugar_carbon_idx:
                lipid_carbon_idx = carbon_idx
                break
              
            is_lipid_carbon = False

            if lipid_carbon_idx != -1:
              if is_fatty_acid_chain(mol, lipid_carbon_idx) or is_glycerol(mol, lipid_carbon_idx) or is_ceramide(mol, lipid_carbon_idx):
                is_lipid_carbon = True


            if is_lipid_carbon:
                glycosidic_matches = True
                break

    if not glycosidic_matches:
        return False, "No glycosidic linkage from sugar ring to a lipid found"

    # 4. Check if molecule contains a phosphate or a sulphate. These will be ignored in the check below.
    has_phosphate = any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms())
    has_sulphate = any(atom.GetAtomicNum() == 16 for atom in mol.GetAtoms() if any(n.GetAtomicNum() == 8 for n in atom.GetNeighbors() ))
    
    # 5. Check for molecular weight and number of atoms (if not phosphate or sulphate)
    if not has_phosphate and not has_sulphate:
        mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
        if mol_wt < 300:
            return False, "Molecular weight too low for glycolipid"
        
        c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        if c_count < 8:
             return False, "Too few carbons for glycolipid"


    return True, "Glycosidic linkage between a carbohydrate and lipid chains is present."
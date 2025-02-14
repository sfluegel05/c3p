"""
Classifies: CHEBI:73754 thiosugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar is a carbohydrate derivative in which one or more of the oxygens
    or hydroxy groups of the parent carbohydrate is replaced by sulfur or -SR,
    where R can be hydrogen or any group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiosugar, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of sulfur atoms
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]
    if not sulfur_atoms:
      return False, "No sulfur atoms found"
    
    # Check for a carbohydrate-like core (ring with multiple hydroxyls or oxygens)
    # A simplified SMARTS pattern to detect a 5 or 6-membered ring with multiple O's or OHs attached to carbons.
    sugar_ring_pattern = Chem.MolFromSmarts("[C;R5,R6]([OX2H,OX2])[C;R5,R6]([OX2H,OX2])[C;R5,R6]([OX2H,OX2])[C;R5,R6]([OX2H,OX2])([C;R5,R6])([C;R5,R6])")
    if not mol.HasSubstructMatch(sugar_ring_pattern):
        return False, "No carbohydrate-like ring structure detected"
    
    
    # Check if any sulfur is directly connected to a sugar ring carbon, or as part of -SR where R is another atom/group.
    sugar_carbon_pattern = Chem.MolFromSmarts("[C;R5,R6]")
    sugar_carbon_atoms = mol.GetSubstructMatches(sugar_carbon_pattern)

    for s_atom in sulfur_atoms:
      for c_match in sugar_carbon_atoms:
          for c_atom_idx in c_match:
              c_atom = mol.GetAtomWithIdx(c_atom_idx)
              if s_atom.GetIdx() in [neighbor.GetIdx() for neighbor in c_atom.GetNeighbors()]:
                  return True, "Sulfur directly connected to a sugar ring carbon"
              for neighbor in s_atom.GetNeighbors():
                #check for S-R, where R is a carbon or hydrogen
                if neighbor.GetAtomicNum() == 1 or neighbor.GetAtomicNum() == 6:
                     return True, "Sulfur connected to sugar ring carbon via a chain"

    # Check for thio-glycosidic bond, a C-S-C between sugar rings
    #First find sugar rings
    sugar_rings_matches = mol.GetSubstructMatches(sugar_ring_pattern)
    if len(sugar_rings_matches) > 1:
        for ring_match_1 in sugar_rings_matches:
            for ring_match_2 in sugar_rings_matches:
                if ring_match_1 != ring_match_2:
                     for s_atom in sulfur_atoms:
                         neighbor_count = 0
                         for neighbor in s_atom.GetNeighbors():
                             for carbon_idx in ring_match_1:
                                if neighbor.GetIdx() == carbon_idx:
                                    neighbor_count +=1
                             for carbon_idx in ring_match_2:
                                if neighbor.GetIdx() == carbon_idx:
                                  neighbor_count +=1
                         if neighbor_count ==2:
                            return True, "Sulfur forms a thio-glycosidic bond between sugar rings"


    # if all checks failed:
    return False, "Sulfur not directly attached to sugar, or as part of a glycosidic bond."
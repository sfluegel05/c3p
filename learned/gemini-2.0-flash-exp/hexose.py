"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is a six-carbon monosaccharide which in its linear form contains either an aldehyde group at position 1 (aldohexose) or a ketone group at position 2 (ketohexose).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for 6 carbons
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    c_count = len(c_atoms)
    if c_count < 6:
        return False, f"Molecule has less than 6 carbons, it has {c_count}."

    # Check for pyranose (6-membered ring with O)
    pyranose_pattern = Chem.MolFromSmarts("O1[CX4][CX4][CX4][CX4][CX4]1")
    pyranose_matches = mol.GetSubstructMatches(pyranose_pattern)
    if pyranose_matches:
        # Check that the ring has exactly 6 carbons
        for match in pyranose_matches:
             ring_atoms = [mol.GetAtomWithIdx(i) for i in match]
             ring_carbon_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 6)
             if ring_carbon_count == 5:
                  continue;
             elif ring_carbon_count != 6:
                 return False, "Pyranose ring does not have exactly 6 carbons."

        return True, "Molecule contains a pyranose (6-membered ring with O)"


    # Check for furanose (5-membered ring with O)
    furanose_pattern = Chem.MolFromSmarts("O1[CX4][CX4][CX4][CX4]1")
    furanose_matches = mol.GetSubstructMatches(furanose_pattern)
    if furanose_matches:
      for match in furanose_matches:
             ring_atoms = [mol.GetAtomWithIdx(i) for i in match]
             ring_carbon_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 6)
             if ring_carbon_count == 4:
                 continue
             elif ring_carbon_count != 5:
                return False, "Furanose ring does not have exactly 5 carbons"
      
      # Check for an exocyclic carbon
      exocyclic_carbon_pattern = Chem.MolFromSmarts("[CX4]~[OX2]1[CX4][CX4][CX4][CX4]1")
      exocyclic_matches = mol.GetSubstructMatches(exocyclic_carbon_pattern)
      if not exocyclic_matches:
        return False, "Furanose ring without an exocyclic carbon"
      return True, "Molecule contains a furanose (5-membered ring with O), and an extra carbon"
      
    # If it's not a ring, then check for linear chain with a carbonyl
    
    # Look for a chain of 6 or more carbons
    carbon_chain_pattern = Chem.MolFromSmarts("[CX4][CX4][CX4][CX4][CX4][CX4]")
    chain_matches = mol.GetSubstructMatches(carbon_chain_pattern)
    
    if chain_matches:
        
        # Check for aldehyde at C1
        aldehyde_pattern = Chem.MolFromSmarts("[CH1](=O)[C]")
        aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)

        if aldehyde_matches:
            for match in aldehyde_matches:
                first_carbon_idx = match[0]
                first_carbon = mol.GetAtomWithIdx(first_carbon_idx)
                # Check if it is a terminal carbon, this is not a perfect solution
                is_terminal = first_carbon.GetDegree() < 3
                if is_terminal:
                     return True, "Molecule is a linear chain with an aldehyde at C1."

        # Check for ketone at C2
        ketone_pattern = Chem.MolFromSmarts("[CX4][C](=O)[CX4]")
        ketone_matches = mol.GetSubstructMatches(ketone_pattern)

        if ketone_matches:
           for match in ketone_matches:
               
            ketone_carbon = mol.GetAtomWithIdx(match[1])
            neighbors = ketone_carbon.GetNeighbors()
            
            c_count = 0
            for neighbor in neighbors:
                if neighbor.GetAtomicNum() == 6:
                    c_count+=1
            if c_count == 2:
                
                # Check if the carbonyl is at C2
                carbon_chain_pattern_c2 = Chem.MolFromSmarts("[CX4]~[CX4](=[OX1])~[CX4]")
                chain_matches_c2 = mol.GetSubstructMatches(carbon_chain_pattern_c2)

                if chain_matches_c2:
                    return True, "Molecule is a linear chain with a ketone at C2"

    
    return False, "Molecule is not a ring or a linear chain with a carbonyl at the 1 or 2 position."
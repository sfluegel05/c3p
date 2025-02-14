"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
"""
Classifies: fatty acid methyl ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    A fatty acid methyl ester is a fatty acid with the carboxyl group esterified with methanol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid methyl ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the methyl ester group pattern (C(=O)OC)
    methyl_ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])OC")
    methyl_ester_matches = mol.GetSubstructMatches(methyl_ester_pattern)
    
    if len(methyl_ester_matches) == 0:
        return False, "No methyl ester group found"
    
    # Check for the methyl group attached to the oxygen on the ester
    has_methyl = False
    for match in methyl_ester_matches:
       ester_o_idx = -1
       for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() == 'O':
                ester_o_idx = idx
                break
       if ester_o_idx != -1:
        for neighbor in mol.GetAtomWithIdx(ester_o_idx).GetNeighbors():
             if neighbor.GetSymbol() == 'C' and neighbor.GetTotalNumHs() == 3:
                has_methyl = True
                break
       if has_methyl:
        break

    if not has_methyl:
       return False, "Methyl group not found on ester oxygen"

    # Function to count consecutive CH2 groups starting from carbonyl carbon
    def count_ch2_chain(mol, start_atom_idx):
        count = 0
        current_atom_idx = start_atom_idx
        visited_atoms = {current_atom_idx}

        while True:
            
            ch2_found = False
            for neighbor in mol.GetAtomWithIdx(current_atom_idx).GetNeighbors():
                if neighbor.GetIdx() not in visited_atoms and neighbor.GetSymbol() == 'C' and neighbor.GetTotalNumHs() == 2:
                     count += 1
                     visited_atoms.add(neighbor.GetIdx())
                     current_atom_idx = neighbor.GetIdx()
                     ch2_found = True
                     break
            if not ch2_found:
                break #no more ch2 in the chain

        return count

    # Check for fatty acid chains (long carbon chains attached to carbonyl of ester)
    fatty_acid_chain_found = False
    for match in methyl_ester_matches:
         carbonyl_idx = -1
         for idx in match:
              atom = mol.GetAtomWithIdx(idx)
              if atom.GetSymbol() == 'C' and atom.GetTotalValence() == 3:
                carbonyl_idx = idx
                break
         if carbonyl_idx != -1:
            ch2_count = count_ch2_chain(mol, carbonyl_idx)
            if ch2_count >= 3: # fatty acid needs a minimum chain length of 3 consecutive CH2 groups 
              fatty_acid_chain_found = True
              break
    
    if not fatty_acid_chain_found:
          return False, "Missing fatty acid chain"


    return True, "Contains a methyl ester group and a fatty acid chain"
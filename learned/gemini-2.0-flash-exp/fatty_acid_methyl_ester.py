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

    # Look for the methyl ester group pattern (C(=O)OC) including methyl hydrogens
    methyl_ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])OC([H])([H])[H]")
    methyl_ester_matches = mol.GetSubstructMatches(methyl_ester_pattern)
    
    if len(methyl_ester_matches) != 1:
        return False, f"Found {len(methyl_ester_matches)} methyl ester groups, need exactly 1"

    
    # Check for fatty acid chains (long carbon chains attached to carbonyl of ester)
    fatty_acid_chain_found = False
    for match in methyl_ester_matches:
         carbonyl_idx = -1
         for idx in match:
              atom = mol.GetAtomWithIdx(idx)
              if atom.GetSymbol() == 'C' and atom.GetTotalValence() == 3: # looking for carbon that is part of ester carbonyl
                carbonyl_idx = idx
                break
         if carbonyl_idx != -1:
             
            alkyl_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
            
            # get all neighbours for carbonyl, not including ester oxygen
            neighbours = [neighbour for neighbour in mol.GetAtomWithIdx(carbonyl_idx).GetNeighbors() if mol.GetAtomWithIdx(neighbour.GetIdx()).GetSymbol() != 'O' ]
           
            for neighbour in neighbours:
                submol = Chem.PathToSubmol(mol, [carbonyl_idx, neighbour.GetIdx()])
                
                if submol.HasSubstructMatch(alkyl_chain_pattern):
                    fatty_acid_chain_found = True
                    break
            
            if fatty_acid_chain_found:
                break
    
    if not fatty_acid_chain_found:
          return False, "Missing fatty acid chain"


    return True, "Contains a methyl ester group and a fatty acid chain"
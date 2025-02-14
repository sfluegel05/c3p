"""
Classifies: CHEBI:61384 sulfolipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    A sulfolipid contains a sulfonic acid residue linked to a lipid via a carbon-sulfur or carbon-oxygen-sulfur bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfolipid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for the sulfonate/sulfonic acid group (-OS(=O)(=O)O-, -S(=O)(=O)OH, -O-S(=O)(=O)-O or -O-S(=O)(=O)=O) connected to a carbon or oxygen
    # More specific SMARTS patterns to match different bonding arrangements
    sulfur_pattern1 = Chem.MolFromSmarts("[CX4]-[S](=[OX1])(=[OX1])-[OX1,OX2]") # -C-S(=O)(=O)-O, -C-S(=O)(=O)-OH
    sulfur_pattern2 = Chem.MolFromSmarts("[CX4]-[OX2]-[S](=[OX1])(=[OX1])-[OX1,OX2]") # -C-O-S(=O)(=O)-O, -C-O-S(=O)(=O)-OH
    
    sulfur_matches1 = mol.GetSubstructMatches(sulfur_pattern1)
    sulfur_matches2 = mol.GetSubstructMatches(sulfur_pattern2)

    if not sulfur_matches1 and not sulfur_matches2:
        return False, "No sulfonate/sulfonic acid group found."
    
    connecting_atoms = []
    if sulfur_matches1:
       for match in sulfur_matches1:
            for idx in match:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 6: # Carbon
                    connecting_atoms.append(atom)
    if sulfur_matches2:
        for match in sulfur_matches2:
            for idx in match:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 8: #Oxygen
                    connecting_atoms.append(atom)
                
    if not connecting_atoms:
        return False, "No carbon or oxygen directly connected to sulfur in the sulfonate/sulfonic acid group."
    
    # 3. Check if the carbon or oxygen is part of a lipid-like moiety: Look for long chain attached to the carbon
    is_lipid = False
    for connecting_atom in connecting_atoms:
         def count_carbons(atom, visited_atoms=None, chain_length=0):
            if visited_atoms is None:
                visited_atoms = set()
            if atom.GetIdx() in visited_atoms:
                return chain_length
            
            visited_atoms.add(atom.GetIdx())
            
            neighbors = [neighbor for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 6 ]
            if not neighbors:
              return chain_length
            
            max_chain_length = chain_length
            for neighbor in neighbors:
              bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
              if bond and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                  max_chain_length = max(max_chain_length, count_carbons(neighbor, visited_atoms, chain_length + 1))
            return max_chain_length
         
         chain_length = count_carbons(connecting_atom)   

         if chain_length >= 8:
            is_lipid = True
            break

    if not is_lipid:
        return False, "Carbon or oxygen connected to sulfur is not part of a lipid-like moiety with at least 8 carbons"

    
    # 4. Check molecular weight as a secondary filter - do not discard molecules if they are above MW 450 as long as they are lipids.  Removed as it caused issues with small sulfolipids.
    #mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    return True, "Contains a sulfonic acid residue linked to a lipid via a carbon-sulfur or carbon-oxygen-sulfur bond."
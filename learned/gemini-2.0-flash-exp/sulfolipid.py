"""
Classifies: CHEBI:61384 sulfolipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    A sulfolipid contains a sulfonic acid residue linked to a lipid via a carbon-sulfur bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfolipid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for the sulfonate/sulfonic acid group (-OS(=O)(=O)O-, -S(=O)(=O)OH, -O-S(=O)(=O)-O or -O-S(=O)(=O)=O) connected to a carbon.
    sulfur_pattern = Chem.MolFromSmarts("[CX4]~[S,O](=[OX1])(=[OX1])[OX1,OX2]")
    sulfur_matches = mol.GetSubstructMatches(sulfur_pattern)
    if not sulfur_matches:
        return False, "No sulfonate/sulfonic acid group found."

    #Get atoms of sulfur or oxygen directly connected to sulfur
    sulfur_atoms = []
    for match in sulfur_matches:
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 16: #Sulfur
                sulfur_atoms.append(atom)
    if not sulfur_atoms:
        return False, "No sulfur found directly connected to a carbon"

    connected_carbons = []
    for sulfur_atom in sulfur_atoms:
        for neighbor in sulfur_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                bond = mol.GetBondBetweenAtoms(sulfur_atom.GetIdx(), neighbor.GetIdx())
                if bond and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                  connected_carbons.append(neighbor)
    
    if not connected_carbons:
        return False, "No carbon directly connected to sulfur via a single bond."
        
    # 3. Check if the carbon is part of a lipid-like moiety: Look for long chain attached to the carbon
    is_lipid = False
    for carbon_atom in connected_carbons:
        #Check connected chain for at least 8 carbons
      
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
          
        chain_length = count_carbons(carbon_atom)
        
        if chain_length >=8:
            is_lipid = True
            break

    if not is_lipid:
      return False, "Carbon connected to sulfur is not part of a lipid-like moiety with at least 8 carbons"

    # 4. Check molecular weight as a secondary filter - do not discard molecules if they are above MW 450 as long as they are lipids
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)

    return True, "Contains a sulfonic acid residue linked to a lipid via a carbon-sulfur or carbon-oxygen-sulfur bond."
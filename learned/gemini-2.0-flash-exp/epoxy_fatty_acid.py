"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: epoxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    An epoxy fatty acid is a fatty acid containing an epoxide ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for epoxide ring pattern (more precise definition of an epoxide)
    epoxide_pattern = Chem.MolFromSmarts("[C]1[O][C]1")
    epoxide_matches = mol.GetSubstructMatches(epoxide_pattern)
    if not epoxide_matches:
        return False, "No epoxide ring found"

    # Look for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found"

    # Check if the epoxide and the acid are linked by a carbon chain
    # We look for a chain of at least 4 carbons connected to the epoxide and the acid group
    # Create SMARTS pattern for a chain linking the epoxide to the acid. We start from the carbons of the epoxide
    for ep_match in epoxide_matches:
      ep_atoms = [mol.GetAtomWithIdx(i) for i in ep_match]
      for acid_match in acid_matches:
        acid_atoms = [mol.GetAtomWithIdx(i) for i in acid_match]
        for ep_atom in ep_atoms:
            for acid_atom in acid_atoms:
                 if ep_atom.GetSymbol() == "C" and acid_atom.GetSymbol() == "C": # both ends must be carbon
                    
                   # Check if there is a path between the epoxide carbon and the acid carbon
                   path = Chem.GetShortestPath(mol, ep_atom.GetIdx(), acid_atom.GetIdx())
                   if path is None:
                        continue #no path between this epoxide/acid pair, try the next
                   
                   #check path length and atom types within the path. Should consist of carbon chains
                   if len(path) >= 5:  # min chain length 4 carbons.  Also must have an O on one end
                      all_carbons = True
                      for atom_idx in path[1:-1]: # do not check starting and final atoms since those are ep/acid
                            if mol.GetAtomWithIdx(atom_idx).GetSymbol() != "C":
                                all_carbons = False
                                break
                      if all_carbons: # check also if it contains double or triple bonds. If it does, its still a fatty acid
                           
                            #check total chain length
                            n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
                            if n_rotatable < 5: # fatty acids should be at least this long
                                return False, "Chain too short for a fatty acid"

                            #Check the number of carbons
                            c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
                            if c_count < 12: # typical fatty acid
                                return False, "Too few carbons for a fatty acid"


                            return True, "Contains an epoxide ring and a fatty acid chain"
    
    # if we reach here, it means no path was found
    return False, "No fatty acid chain linking the epoxide and the carboxylic acid group"
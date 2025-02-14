"""
Classifies: CHEBI:33567 catecholamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine is a benzene ring with hydroxyl groups at positions 1 and 2 (a catechol),
    and a 2-aminoethyl group attached at position 4, with possible substitutions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catecholamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for catechol moiety (1,2-dihydroxybenzene)
    catechol_pattern = Chem.MolFromSmarts("c1c(O)c(O)c(*)c(*)c1")
    catechol_matches = mol.GetSubstructMatches(catechol_pattern)
    if not catechol_matches:
         return False, "No catechol moiety (1,2-dihydroxybenzene) found"

    # Check for 2-aminoethyl group (-CH2-CH2-NH2)
    aminoethyl_pattern = Chem.MolFromSmarts("[#6]~[#6]~[#7]")
    aminoethyl_matches = mol.GetSubstructMatches(aminoethyl_pattern)
    if not aminoethyl_matches:
        return False, "No 2-aminoethyl group found"
    
    # Check that the nitrogen is not directly attached to the 1 or 2 carbons of the benzene ring of the catechol moiety
    nitrogen_attachment_pattern = Chem.MolFromSmarts("[#7]~[#6][c;!$(c[O]);!$(c[O])]")
    if not mol.HasSubstructMatch(nitrogen_attachment_pattern):
      return False, "2-aminoethyl group not correctly attached to the benzene ring"
    
    # Check that at least one of the aminoethyl groups is attached to the catechol ring
    found_valid_attachment = False
    for aminoethyl_match in aminoethyl_matches:
        nitrogen_atom = mol.GetAtomWithIdx(aminoethyl_match[2]) #get the nitrogen atom
        
        for catechol_match in catechol_matches:
            
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in catechol_match]
            
            #Check if any of the two carbons attached to the nitrogen is a ring atom that IS NOT the hydroxyl attached atoms.
            
            for neighbor in nitrogen_atom.GetNeighbors():
                if neighbor.GetIdx() != aminoethyl_match[1]: #get atom that attaches to nitrogen
                    if neighbor in ring_atoms:
                        
                        is_attached_to_correct_position = True
                        #check that the nitrogen is not directly connected to an oxygenated carbon of the catechol
                        for ring_atom in ring_atoms:
                            if (ring_atom.GetIdx() == catechol_match[1] or ring_atom.GetIdx() == catechol_match[2]): # these are the carbon attacehd to the hydroxyls, we do not want direct bond
                                for neighbor_of_nitrogen in nitrogen_atom.GetNeighbors():
                                    if neighbor_of_nitrogen == ring_atom:
                                        is_attached_to_correct_position = False
                                        break
                            if not is_attached_to_correct_position:
                                break
                        if is_attached_to_correct_position:
                            found_valid_attachment = True
                            break
            if found_valid_attachment:
              break
        if found_valid_attachment:
          break

    if not found_valid_attachment:
        return False, "2-aminoethyl group not correctly attached to the benzene ring"
    

    return True, "Contains a catechol moiety and a 2-aminoethyl group attached to the ring"
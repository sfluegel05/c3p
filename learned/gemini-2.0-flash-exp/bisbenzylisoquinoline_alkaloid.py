"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
"""
Classifies: bisbenzylisoquinoline alkaloid
"""
from rdkit import Chem
from rdkit.Chem import AllChem


def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    A bisbenzylisoquinoline alkaloid is characterized by two benzylisoquinoline units linked by ether bridges.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
         bool: True if molecule is a bisbenzylisoquinoline alkaloid, False otherwise
         str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core benzylisoquinoline substructure. The nitrogen can be charged.
    # The [CX4] is benzylic carbon, next to the isoquinoline ring.
    # The [!H] means 'non-hydrogen atom', so it can be any atom connected to the benzylic C
    benzylisoquinoline_pattern = Chem.MolFromSmarts("[NX3,NX4+]1([CH2][!H])=C2C=CC=CC=2c3ccccc3")

    if benzylisoquinoline_pattern is None:
        return None, "Invalid SMARTS Pattern"

    # Find all matches of the benzylisoquinoline pattern
    matches = mol.GetSubstructMatches(benzylisoquinoline_pattern)
    if len(matches) < 2:
         return False, f"Found {len(matches)} benzylisoquinoline units, need at least 2"


    # Check for ether linkages (C-O-C)
    ether_pattern = Chem.MolFromSmarts("[CX4][OX2][CX4]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)


    # If no ether linkages, it might still be a bisbenzylisoquinoline, but with another connection,
    # such as a carbon-carbon bridge, a methylenedioxy bridge or a direct bond.
    # The connection can only occur between two benzylisoquinoline fragments.
    if not ether_matches:
        methylenedioxy_pattern = Chem.MolFromSmarts("C[O][C][O]C")
        methylenedioxy_matches = mol.GetSubstructMatches(methylenedioxy_pattern)

        if not methylenedioxy_matches:
           
            # Check for direct carbon-carbon bonds connecting two benzylisoquinolines.
            # the following checks if two benzylisoquinoline fragments are connected
            # by any bond.
            is_directly_connected = False
            for match1 in matches:
                 for match2 in matches:
                     if match1 == match2:
                         continue
                     for idx1 in match1:
                         atom1 = mol.GetAtomWithIdx(idx1)
                         for idx2 in match2:
                            atom2 = mol.GetAtomWithIdx(idx2)
                            if mol.GetBondBetweenAtoms(atom1.GetIdx(), atom2.GetIdx()):
                                 is_directly_connected = True
                                 break
                         if is_directly_connected:
                             break
                     if is_directly_connected:
                         break


            if not is_directly_connected:
                return False, "No ether, methylenedioxy or direct bond connecting two benzylisoquinoline units"

        else:
            # We found a methylenedioxy bridge, now we check if this bridges two bzisoquinolines
            is_methylenedioxy_connected = False
            for bridge_match in methylenedioxy_matches:
                bridge_atoms = [mol.GetAtomWithIdx(idx) for idx in bridge_match]
                for bzisoq_match in matches:
                     
                    for idx in bzisoq_match:
                        bzisoq_atom = mol.GetAtomWithIdx(idx)
                        
                        for bridge_atom in bridge_atoms:
                            if bridge_atom.GetIdx() in [neigh.GetIdx() for neigh in bzisoq_atom.GetNeighbors()]:
                                 for bzisoq_match2 in matches:
                                     if bzisoq_match == bzisoq_match2:
                                        continue
                                     for idx2 in bzisoq_match2:
                                         bzisoq_atom2 = mol.GetAtomWithIdx(idx2)
                                         for bridge_atom2 in bridge_atoms:
                                              if bridge_atom2.GetIdx() in [neigh.GetIdx() for neigh in bzisoq_atom2.GetNeighbors()]:
                                                   is_methylenedioxy_connected = True
                                                   break
                                     if is_methylenedioxy_connected:
                                        break
                        if is_methylenedioxy_connected:
                            break
                    if is_methylenedioxy_connected:
                        break
            if not is_methylenedioxy_connected:
                 return False, "The methylenedioxy bridge does not connect the two benzylisoquinoline units"
            
        return True, "Contains two benzylisoquinoline units connected by at least a methylenedioxy bridge or a direct bond"



    #Check that the ether is between two benzylisoquinolines
    ether_count = 0
    for match in ether_matches:
         atom1 = mol.GetAtomWithIdx(match[0])
         atom2 = mol.GetAtomWithIdx(match[1])
         atom3 = mol.GetAtomWithIdx(match[2])
         
         is_connected_to_bzisoquinolines = False
         for bzisoq_match in matches:
            for bzisoq_idx in bzisoq_match:
               bzisoq_atom = mol.GetAtomWithIdx(bzisoq_idx)
               if atom1.GetIdx() in [neighbor.GetIdx() for neighbor in bzisoq_atom.GetNeighbors()]:
                   
                    for bzisoq_match2 in matches:
                        if bzisoq_match == bzisoq_match2:
                            continue
                        for bzisoq_idx2 in bzisoq_match2:
                             bzisoq_atom2 = mol.GetAtomWithIdx(bzisoq_idx2)
                             if atom3.GetIdx() in [neighbor.GetIdx() for neighbor in bzisoq_atom2.GetNeighbors()]:
                                 is_connected_to_bzisoquinolines = True
                                 break
                        if is_connected_to_bzisoquinolines:
                            break
            if is_connected_to_bzisoquinolines:
                ether_count+=1
         
    if ether_count < 1:
         return False, "The ether bridge does not connect the two benzylisoquinoline units"


    return True, "Contains two benzylisoquinoline units connected by at least one ether bridge"
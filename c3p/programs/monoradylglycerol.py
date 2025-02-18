"""
Classifies: CHEBI:76575 monoradylglycerol
"""
"""
Classifies: monoradylglycerol
Definition: Any lipid that is glycerol bearing a single acyl, alkyl or alk-1-enyl substituent at an unspecified position.
In our algorithm, we search for a glycerol backbone (a series of three contiguous sp3 carbons) in which 
each carbon bears one oxygen substituent. In a free glycerol all three oxygens would be –OH; here we demand 
that exactly one oxygen lacks a hydrogen (“substituted” by an acyl/alkyl group), while the other two are free –OH.
Furthermore, the substituent connected via the substituted oxygen is required to be a carbon-based chain of minimum length.
"""

from rdkit import Chem

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    A monoradylglycerol is defined as glycerol bearing a single acyl/alkyl/alk-1-enyl substituent.
    In our model, we look for a three-carbon glycerol-like backbone where each carbon bears one oxygen.
    Exactly one of these oxygens must be substituted (i.e. not free OH) and its attached group must be carbon-based 
    with a minimum chain length (here >= 4 carbons).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a monoradylglycerol, False otherwise.
        str: A reason for the classification.
    """
    
    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Helper function: recursively compute the longest chain of carbons starting from a given carbon atom.
    def longest_carbon_chain(atom, parent_idx, visited):
        # Only count carbon atoms
        if atom.GetAtomicNum() != 6:
            return 0
        visited.add(atom.GetIdx())
        max_length = 1  # count current atom
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() == parent_idx:
                continue
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                length = 1 + longest_carbon_chain(nbr, atom.GetIdx(), visited.copy())
                if length > max_length:
                    max_length = length
        return max_length

    # Iterate over atoms to find a candidate for the central (middle) carbon of the glycerol backbone.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        # Only consider sp3 carbons (typical for glycerol backbone)
        if atom.GetHybridization().name != "SP3":
            continue
        
        # Identify sp3 carbon neighbors of the candidate atom.
        carbon_neighbors = [nbr for nbr in atom.GetNeighbors() 
                            if nbr.GetAtomicNum() == 6 and nbr.GetHybridization().name == "SP3"]
        if len(carbon_neighbors) < 2:
            continue
        
        # Consider each pair of neighboring carbons as potential backbone ends (c1 and c3).
        for i in range(len(carbon_neighbors)):
            for j in range(i+1, len(carbon_neighbors)):
                c1 = carbon_neighbors[i]
                c3 = carbon_neighbors[j]
                # For a linear glycerol backbone, we expect that c1 and c3 are not directly bonded.
                if mol.GetBondBetweenAtoms(c1.GetIdx(), c3.GetIdx()):
                    continue  # skip if they are directly bonded (which might indicate a ring or a branched system)
                
                # Candidate backbone: c1 - atom (middle) - c3
                backbone_atoms = [c1, atom, c3]
                # For each backbone carbon, look for exactly one oxygen neighbor that is not part of the backbone.
                oxygen_neighbors = []
                valid_backbone = True
                for carbon in backbone_atoms:
                    # Gather oxygen neighbors not in the backbone.
                    o_neighbors = [nbr for nbr in carbon.GetNeighbors() 
                                   if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in [a.GetIdx() for a in backbone_atoms]]
                    if len(o_neighbors) != 1:
                        valid_backbone = False
                        break
                    oxygen_neighbors.append(o_neighbors[0])
                if not valid_backbone:
                    continue
                
                # Count free hydroxyl oxygens and those that are substituted.
                # We assume a free -OH will have at least one hydrogen (implicit H counts via GetTotalNumHs).
                free_OH_count = 0
                substituted_index = -1
                for idx, o_atom in enumerate(oxygen_neighbors):
                    if o_atom.GetTotalNumHs() > 0:
                        free_OH_count += 1
                    else:
                        substituted_index = idx
                # We require exactly two free OH groups and one substituted oxygen.
                if free_OH_count != 2 or substituted_index < 0:
                    continue
                
                # Now examine the substituent group attached via the substituted oxygen.
                sub_oxygen = oxygen_neighbors[substituted_index]
                # For this oxygen, choose neighbors that are not the associated backbone carbon.
                backbone_ids = [a.GetIdx() for a in backbone_atoms]
                sub_neighbors = [nbr for nbr in sub_oxygen.GetNeighbors() if nbr.GetIdx() not in backbone_ids]
                if not sub_neighbors:
                    continue
                # Ideally the substituent is carbon-based.
                sub_atom = sub_neighbors[0]
                if sub_atom.GetAtomicNum() != 6:
                    continue
                # Calculate the longest carbon chain starting from the substituent atom.
                chain_length = longest_carbon_chain(sub_atom, sub_oxygen.GetIdx(), set())
                min_chain_length = 4  # heuristic minimum chain length for a valid substituent
                if chain_length < min_chain_length:
                    continue
                
                reason = (f"Found glycerol backbone (atoms {', '.join(str(a.GetIdx()) for a in backbone_atoms)}) with substitution "
                          f"at position {substituted_index+1}; substituent chain length = {chain_length}, free -OH groups = 2.")
                return True, reason

    return False, "No glycerol backbone with a single substituted oxygen (and two free OH groups) found."

# Example usage and testing of several SMILES strings representing monoradylglycerols:
if __name__ == '__main__':
    test_smiles = [
        "O=C(OC[C@@H](O)CO)CCCCCCC/C=C(\\CCCCCCCC)/C",  # 2,3-dihydroxypropyl (Z)-10-methyloctadec-9-enoate
        "CCCCCCCC(=O)OCC(O)CO",                         # 1-monooctanoylglycerol
        "CCCCCCCCCCCCCCCCCC(=O)OC[C@H](O)CO",           # 3-stearoyl-sn-glycerol
        "O(C(=O)CCCCCCCCCCC/C=C\\C/C=C\\CCCCC)C[C@@H](O)CO",  # MG(22:2(13Z,16Z)/0:0/0:0)
        "O(C(=O)CCCCCCCCCCCCCCCCCCCCC)C[C@@H](O)CO",     # MG(0:0/22:0/0:0)
        "CCCCCCCCCCCCCCCCCCCOC[C@@H](O)CO",              # 1-O-octadecyl-sn-glycerol
    ]
    for s in test_smiles:
        result, reason = is_monoradylglycerol(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")
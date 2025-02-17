"""
Classifies: CHEBI:26167 polar amino acid
"""
#!/usr/bin/env python
"""
Classifies: Polar Amino Acid
Definition: Any amino acid whose side chain is capable of forming one or more hydrogen bonds.
This implementation identifies an amino acid backbone via SMARTS and then extracts
the side chain for a search for polar functional groups (oxygen, nitrogen, or –SH).
Note: Glycine, which has only a hydrogen as a side chain, is not considered polar.
"""

from rdkit import Chem

def is_polar_amino_acid(smiles: str):
    """
    Determines whether the input molecule (SMILES string) is a polar amino acid.
    
    The procedure is as follows:
    1. Use a SMARTS pattern to detect the amino acid backbone.
       The SMARTS used is: "[N:1][C:2]([*:3])C(=O)[O:4]"
         - [N:1] is the amino group.
         - [C:2] is the α‐carbon.
         - [*:3] is the side chain attachment (which could be any atom).
         - C(=O)[O:4] is the carboxyl group (note the extra unlabeled C in the pattern).
    2. If the amino acid backbone is detected and there is only one match (to avoid peptides),
       extract the indices corresponding to the labeled atoms.
    3. Do a graph traversal starting from the side chain atom (labeled as 3) but do not go
       back to the α‐carbon. Then, search among the side chain atoms for polar atoms: oxygen,
       nitrogen, or a sulfur carrying at least one hydrogen (–SH).
    4. Return True (and a reason) if a polar group is found; otherwise, return False.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): True with an explanation if the molecule is a polar amino acid,
                     False with a reason if it is not.
    """
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS to capture the amino acid backbone.
    # The pattern has an unlabeled C in the carboxyl group; hence overall 5 atoms.
    aa_smarts = "[N:1][C:2]([*:3])C(=O)[O:4]"
    patt = Chem.MolFromSmarts(aa_smarts)
    # Get all matching substructures in the molecule
    matches = mol.GetSubstructMatches(patt)
    
    if not matches:
        return False, "Amino acid backbone (N-αC-C(=O)O) not detected"
    
    # If more than one amino acid backbone is found, assume the molecule is a peptide or larger.
    if len(matches) > 1:
        return False, "Multiple amino acid backbones detected; molecule appears to be a peptide or larger compound"
    
    # Only one match found; however, the match tuple length is 5 (including the unlabeled atom).
    match = matches[0]
    
    # Create a mapping: for every query atom in the SMARTS pattern that has a molAtomMapNumber,
    # map that number (as a string) to the corresponding atom index in the molecule.
    labeled_mapping = {}
    # Note: patt.GetAtoms() gives the atoms in the query in order.
    for i, atom in enumerate(patt.GetAtoms()):
        if atom.HasProp("molAtomMapNumber"):
            mapnum = atom.GetProp("molAtomMapNumber")
            labeled_mapping[mapnum] = match[i]
    
    # Ensure that the mapping has all 4 required entries.
    if not all(x in labeled_mapping for x in ["1", "2", "3", "4"]):
        return False, "SMARTS mapping incomplete, backbone not properly detected"
    
    idx_n = labeled_mapping["1"]      # Nitrogen of amino group
    idx_alpha = labeled_mapping["2"]  # α‐Carbon
    idx_side = labeled_mapping["3"]   # Side chain branch atom
    # idx_o = labeled_mapping["4"]     # Carboxyl oxygen (not used further)
    
    # Check for glycine: if the side chain is just a hydrogen, it's typically considered non‐polar.
    side_atom = mol.GetAtomWithIdx(idx_side)
    if side_atom.GetAtomicNum() < 6:
        return False, "Side chain is hydrogen (glycine), thus not polar by our definition"
    
    # Perform a depth-first traversal of the side chain subgraph.
    # Start at the side chain branch and do not traverse back into the backbone (i.e. the α‐carbon).
    sidechain_atoms = set()
    stack = [idx_side]
    while stack:
        cur_idx = stack.pop()
        if cur_idx in sidechain_atoms:
            continue
        sidechain_atoms.add(cur_idx)
        cur_atom = mol.GetAtomWithIdx(cur_idx)
        for nb in cur_atom.GetNeighbors():
            if nb.GetIdx() == idx_alpha:
                continue  # do not go back to the backbone
            if nb.GetIdx() not in sidechain_atoms:
                stack.append(nb.GetIdx())
                
    # Analyze the collected side chain atoms for polar elements or polar groups.
    polar_found = False
    polar_features = []
    for aidx in sidechain_atoms:
        atom = mol.GetAtomWithIdx(aidx)
        anum = atom.GetAtomicNum()
        # Oxygen (atomic number 8) or Nitrogen (atomic number 7) are polar.
        if anum == 8 or anum == 7:
            polar_found = True
            polar_features.append(atom.GetSymbol())
        elif anum == 16:
            # For sulfur, we require at least one hydrogen attached (–SH) to be polar.
            if atom.GetTotalNumHs() > 0:
                polar_found = True
                polar_features.append("SH")
    
    if polar_found:
        return True, f"Side chain contains polar feature(s): {', '.join(polar_features)}"
    else:
        return False, "Side chain does not contain a polar functional group capable of hydrogen bonding"

# Example usage (for testing purposes):
if __name__ == "__main__":
    test_smiles = [
        "NC(CCC(N)=O)C(O)=O",        # glutamine, polar
        "NC(CO)C(O)=O",             # serine, polar
        "N[C@H](CS)C(O)=O",          # D-cysteine, polar due to –SH
        "N[C@@H](Cc1c[nH]cn1)C(O)=O", # L-histidine, polar (N in side chain)
        "N[C@@H](CCC(N)=O)C(O)=O",    # L-glutamine, polar
        "NCCCC[C@@H](N)C(O)=O"        # D-lysine, polar (side chain N)
    ]
    for s in test_smiles:
        result, reason = is_polar_amino_acid(s)
        print(f"SMILES: {s}\nIs polar amino acid? {result} ({reason})\n")
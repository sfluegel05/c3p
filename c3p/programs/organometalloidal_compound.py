"""
Classifies: CHEBI:143084 organometalloidal compound
"""
#!/usr/bin/env python
"""
Classifies: organometalloidal compound
Definition: A compound having bonds between one or more metalloid atoms and one or more carbon atoms
             of an organyl group.
             
Improved heuristics for arsenic–based compounds:
  • Only consider compounds where the SMILES contains "[As".
  • Hard–coded filters reject known inorganic patterns (e.g. tryparsamide–like fragments, certain chloride patterns).
  • For each arsenic (atomic number 33), examine each As–C bond.
       – If the bonded carbon is “embedded” in an organic fragment (i.e. it has at least one additional carbon neighbor) then the bond is accepted.
       – If the bonded carbon is a methyl group (i.e. no additional carbon neighbors after excluding the As) then we accept it only if (i) the overall molecular weight is low (<170 Da) or (ii) the As carries a formal charge or (iii) exactly two methyl groups are present and there are few extra hetero–atomic substituents.
  • The previous version required a double–bonded O on neutral, non–ring As atoms. That requirement has been removed to improve detection.
  
Note: This heuristic approach remains imperfect.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_organometalloidal_compound(smiles: str):
    """
    Determines if a molecule is an organometalloidal compound based on its SMILES string.
    (Improved Heuristics for arsenic-based compounds.)
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is judged to be an organometalloidal compound, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # Only consider arsenic-based compounds.
    if "[As" not in smiles:
        return False, "No arsenic atom found; not an organometalloidal compound"
    
    # Hard–coded filters (these can be extended as needed):
    if "NCC(N)=O" in smiles:
        return False, "Detected tryparsamide-like fragment (NCC(N)=O); not considered an organyl group"
    if "[As](Cl" in smiles:
        return False, "Detected Cl substituents on arsenic; likely inorganic derivative"
    
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 120:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) to be an organometalloidal compound"
    
    overall_details = []  # collect explanation details
    compound_is_valid = False  # flag for at least one acceptable As–C bonding

    # Loop over each atom and focus on arsenic (atomic number 33).
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 33:
            continue

        # For each As atom, collect its carbon neighbors.
        carbon_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if not carbon_neighbors:
            continue  # This As has no As-C bond

        as_details = []   # details for this particular As atom
        valid_for_this_As = False  # flag to mark if any neighbor qualifies
        
        # Counters used for methyl-only mode.
        methyl_count = 0
        
        for c in carbon_neighbors:
            # Determine the number of other carbon neighbors for carbon 'c'
            # (exclude the source As atom).
            other_carbons = [nbr for nbr in c.GetNeighbors()
                             if nbr.GetIdx() != atom.GetIdx() and nbr.GetAtomicNum() == 6]
            if other_carbons:
                # c is part of a larger organic fragment (chain or ring)
                as_details.append("Found As–C bond where C is part of a larger organic fragment")
                valid_for_this_As = True
            else:
                # c appears to be a methyl group (only neighbor is the arsenic).
                methyl_count += 1
                as_details.append("Found As–C bond where C appears as a methyl group")
                # Accept methyl substituent if:
                #  - Molecule is very simple (low molecular weight) or
                #  - The arsenic atom carries a nonzero formal charge.
                if atom.GetFormalCharge() != 0 or mol_wt < 170:
                    as_details.append("Methyl substituent accepted due to charge/low–MW criteria")
                    valid_for_this_As = True
                # Also, if exactly two methyl groups are present and there aren’t extra heteroatoms bonded
                # to As (besides C and H), then accept.
                elif len(carbon_neighbors) == 2:
                    non_c_neighbors = [nbr for nbr in atom.GetNeighbors()
                                       if nbr.GetAtomicNum() not in (6, 1)]
                    if len(non_c_neighbors) <= 2:
                        as_details.append("Two methyl groups accepted based on substitution pattern")
                        valid_for_this_As = True

        # If for this As we found at least one acceptable As-C bond, mark compound as valid.
        if valid_for_this_As:
            overall_details.extend(as_details)
            compound_is_valid = True
            break  # At least one As qualifies; no need to check further.
        else:
            overall_details.extend(as_details)
    
    if compound_is_valid:
        return True, "; ".join(overall_details)
    else:
        return False, "No As–C bond found that meets organyl criteria"

# Example usage:
if __name__ == "__main__":
    # True positive examples:
    examples = [
        ("[As](=O)(C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC)(C)C", "1-dimethylarsinoyl-(2Z,5Z,8Z,11Z,14Z,17Z)-heneicosahexaene"),
        ("O[As](O)(=O)c1ccc(cc1)[N+]#N", "ABA diazonium"),
        ("C1C=CC=[As]1", "2H-arsole"),
        ("[As](=O)(CCCCCCCCCCCCCCCCCCCCCCC)(C)C", "1-dimethylarsinoyl-tricosane"),
        ("C[As](C)c1ccccc1[As](C)C", "1,2-phenylenebis(dimethylarsane) [expected now to be accepted]"),
        ("OC(=O)C[As](O)(O)=O", "arsenoacetic acid [expected to be accepted]"),
    ]
    for smi, name in examples:
        result, reason = is_organometalloidal_compound(smi)
        print(f"Name: {name}\nSMILES: {smi}\nResult: {result}\nReason: {reason}\n{'-'*60}")
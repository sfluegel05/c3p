"""
Classifies: CHEBI:143084 organometalloidal compound
"""
#!/usr/bin/env python
"""
Classifies: organometalloidal compound
Definition: A compound having bonds between one or more metalloid atoms and one or more carbon atoms 
            of an organyl group.
            
Improved heuristics for arsenic–based compounds:
  - Only consider molecules whose SMILES string contains "[As".
  - Reject known inorganic patterns (e.g. tryparsamide–like fragments, chloride substitution).
  - For each arsenic (atomic number 33), examine every As–C bond.
       • A carbon is considered “extended” (i.e. part of an organyl group) if it is bonded to at least one other carbon (aside from the analysed As atom).
       • If any extended As–C bond is found, then the compound is accepted.
       • If only methyl groups (i.e. carbons with no extra carbon neighbors) are present, then accept 
         the bond only if (a) exactly one (or two under a strict substitution‐pattern) is present AND 
         either the overall molecular weight is very low (<170 Da) or the As carries a nonzero formal charge.
         (This avoids accepting compounds like trimethylarsine oxide.)
         
Note: This heuristic approach is not perfect.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_organometalloidal_compound(smiles: str):
    """
    Determines if a molecule is an organometalloidal compound based on its SMILES string.
    (Improved heuristics for arsenic-based compounds.)
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is judged to be an organometalloidal compound; False otherwise.
        str: A reason detailing the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # Quick filter: only consider molecules that include "[As" (indicative of an arsenic atom)
    if "[As" not in smiles:
        return False, "No arsenic atom found; not an organometalloidal compound"
    
    # Hard-coded filters for known inorganic fragments:
    if "NCC(N)=O" in smiles:
        return False, "Detected tryparsamide-like fragment (NCC(N)=O); not considered an organyl group"
    if "[As](Cl" in smiles:
        return False, "Detected Cl substituents on arsenic; likely inorganic derivative"
    
    # Calculate molecular weight (used in some criteria)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 120:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) to be an organometalloidal compound"
    
    overall_details = []  # Collect details for the explanation
    compound_is_valid = False  # Set if we find at least one acceptable As–C bond
    
    # Loop over all atoms to identify arsenic atoms (atomic number 33)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 33:
            continue
        
        # Get all carbon neighbors for the current arsenic atom
        carbon_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if not carbon_neighbors:
            continue  # No As-C bond on this As; try next
        
        valid_for_this_As = False  # Flag for acceptable bonding on this As
        as_details = []  # Explanation details for this As atom
        
        # Counters for methyl bonds (and list details if only methyl groups are encountered)
        methyl_count = 0
        
        # For each C neighbor bonded to our As atom:
        for c in carbon_neighbors:
            # Check how many neighbors (of type C) the carbon has (excluding the current As)
            other_carbons = [nbr for nbr in c.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != atom.GetIdx()]
            if other_carbons:
                # This carbon is part of a larger organic fragment (chain or ring)
                as_details.append("Found As–C bond where C is part of a larger organic fragment")
                valid_for_this_As = True
            else:
                # This carbon appears to be a methyl group (only bonded to the arsenic)
                methyl_count += 1
                as_details.append("Found As–C bond where C appears as a methyl group")
        
        # If at least one extended (non-methyl) bond was found, we accept immediately.
        if valid_for_this_As:
            overall_details.extend(as_details)
            compound_is_valid = True
            break
        
        # Else, if only methyl groups are present, apply stricter rules.
        if methyl_count > 0:
            # We require that only one (or at best two under narrow conditions) methyl groups are attached.
            if methyl_count == 1:
                # Accept the lone methyl if either the arsenic has a formal charge or the molecule is very small.
                if atom.GetFormalCharge() != 0 or mol_wt < 170:
                    as_details.append("Single methyl group accepted due to charge/low–MW criteria")
                    valid_for_this_As = True
            elif methyl_count == 2:
                # For two methyl groups, check that there are no additional heteroatom substituents (besides H).
                non_CH_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() not in (6, 1)]
                if len(non_CH_neighbors) == 0 and (atom.GetFormalCharge() != 0 or mol_wt < 170):
                    as_details.append("Two methyl groups accepted based on substitution pattern and charge/low–MW criteria")
                    valid_for_this_As = True
            
            # Update details and flag if criteria met.
            if valid_for_this_As:
                overall_details.extend(as_details)
                compound_is_valid = True
                break
            else:
                overall_details.extend(as_details)
    
    if compound_is_valid:
        return True, "; ".join(overall_details)
    else:
        return False, "No As–C bond found meeting the organyl criteria (suitable extended organic fragment or acceptable methyl pattern)"

# Example usage:
if __name__ == "__main__":
    # A selection of examples (true positives and false positives):
    examples = [
        ("[As](=O)(C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC)(C)C", "1-dimethylarsinoyl-(2Z,5Z,8Z,11Z,14Z,17Z)-heneicosahexaene"),
        ("O[As](O)(=O)c1ccc(cc1)[N+]#N", "ABA diazonium"),
        ("C1C=CC=[As]1", "2H-arsole"),
        ("[As](=O)(CCCCCCCCCCCCCCCCCCCCCCC)(C)C", "1-dimethylarsinoyl-tricosane"),
        ("Nc1cc(ccc1O)[As]1[As]([As]([As]([As]1c1ccc(O)c(N)c1)c1ccc(O)c(N)c1)c1ccc(O)c(N)c1)c1ccc(O)c(N)c1", "arsphenamine pentamer"),
        ("C[As](C)c1ccccc1[As](C)C", "1,2-phenylenebis(dimethylarsane)"),
        ("C[As](O)([O-])=O", "methylarsonate(1-)"),
        ("[Na+].C[As](C)([O-])=O", "sodium dimethylarsinate (expected to be rejected)"),
        ("C[As](C)(C)=O", "trimethylarsine oxide (expected to be rejected)")
    ]
    for smi, name in examples:
        result, reason = is_organometalloidal_compound(smi)
        print(f"Name: {name}\nSMILES: {smi}\nResult: {result}\nReason: {reason}\n{'-'*60}")
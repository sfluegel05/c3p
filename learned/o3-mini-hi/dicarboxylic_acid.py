"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
"""
Classifies: Any carboxylic acid containing two carboxy groups (dicarboxylic acid)

Improvement in this version:
  • Instead of relying solely on SMARTS substructure matches, we loop over all carbon atoms to 
    “manually” detect a carboxylic acid moiety: a carbon bound to one oxygen by a double bond 
    and one oxygen by a single (acidic) bond.
  • We also check that the candidate acid carbon has no nitrogen neighbors (to avoid, for example, 
    matching carboxamide groups in peptides).
  • Finally, since many peptides are being falsely classified as dicarboxylic acids (even if they have 
    two terminal acid groups) we reject any molecule containing an amide bond (C(=O)N).
    
A dicarboxylic acid is accepted only when exactly two such (free) acid groups are present.
"""

from rdkit import Chem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid (exactly two free carboxy acid groups)
    based on its SMILES string.
    
    The algorithm:
      1. Parses the SMILES string and adds explicit hydrogens.
      2. Checks for the presence of amide bonds (C(=O)N); if any are found, the molecule is assumed
         to be a peptide or similar, and is not classified as a simple dicarboxylic acid.
      3. Iterates over all carbon atoms and for each checks for:
           - Exactly one double bond to an oxygen (i.e. a C=O) AND
           - Exactly one single bond to an oxygen that is “acidic” (i.e. it either carries at least one 
             hydrogen or bears a –1 formal charge).
           - The candidate carbon must not be directly bonded to any nitrogen.
      4. If exactly two unique acid carbons are found, the molecule is classified as a dicarboxylic acid.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a dicarboxylic acid, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that acidic -OH groups are visible.
    mol = Chem.AddHs(mol)
    
    # Reject molecules that contain an amide bond (C(=O)N), which often signal peptides.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if amide_pattern is not None and mol.HasSubstructMatch(amide_pattern):
        return False, "Molecule contains amide bonds, likely a peptide rather than a dicarboxylic acid"
    
    acid_carbons = set()  # to record indices of carbons that are part of a free carboxy group

    # Iterate over all atoms; only carbon atoms can be the central acid carbon.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue  # not a carbon
        
        # Check neighbors to see if we have exactly:
        #   - One oxygen attached with a double bond (C=O)
        #   - One oxygen attached with a single bond that is acidic (has at least one hydrogen or is negatively charged)
        double_count = 0
        single_acidic = 0
        # Also, if the carbon is directly bonded to any nitrogen, we want to skip it.
        has_nitrogen = False
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 7:
                has_nitrogen = True
                break
        if has_nitrogen:
            continue

        # Loop over bonds from the carbon atom and check bonds to oxygens
        for bond in atom.GetBonds():
            # Get the neighboring atom (other than the current carbon)
            o_atom = bond.GetOtherAtom(atom)
            if o_atom.GetAtomicNum() != 8:
                continue
            # Check bond type
            btype = bond.GetBondType()
            if btype == Chem.BondType.DOUBLE:
                double_count += 1
            elif btype == Chem.BondType.SINGLE:
                # Consider the oxygen acidic if it has at least one hydrogen explicitly,
                # or if its formal charge is -1.
                if o_atom.GetTotalNumHs() > 0 or o_atom.GetFormalCharge() == -1:
                    single_acidic += 1

        # In a proper carboxylic acid group, we expect one double and one “acidic” single bond.
        if double_count == 1 and single_acidic == 1:
            acid_carbons.add(atom.GetIdx())

    num_acid_groups = len(acid_carbons)
    if num_acid_groups == 2:
        return True, "Molecule contains exactly two carboxylic acid groups"
    elif num_acid_groups < 2:
        return False, f"Found only {num_acid_groups} carboxylic acid group(s); need exactly two."
    else:
        return False, f"Found {num_acid_groups} carboxylic acid groups; dicarboxylic acid must have exactly two."

# Example test cases for manual checking (should be run under __main__)
if __name__ == "__main__":
    examples = {
        "fumarylacetic acid": "OC(=O)CC(=O)\\C=C\\C(O)=O",
        "nocardicin A": "N[C@H](CCOc1ccc(cc1)C(=N\\O)\\C(=O)N[C@H]1CN([C@@H](C(O)=O)c2ccc(O)cc2)C1=O)C(O)=O",
        "meso-tartaric acid": "O[C@@H]([C@@H](O)C(O)=O)C(O)=O",
        "3-ethylmalic acid": "CCC(C(O)C(O)=O)C(O)=O",
        "(4Z)-2-oxohept-4-enedioic acid": "C(=C\\CC(C(O)=O)=O)\\CC(O)=O",
        "Ile-Ser-Asp (peptide false positive?)": "O=C(N[C@H](C(=O)N[C@@H](CC(O)=O)C(O)=O)CO)[C@@H](N)[C@H](CC)C",
        "2-aminoadipic acid": "NC(CCCC(O)=O)C(O)=O",
    }
    for name, smi in examples.items():
        res, reason = is_dicarboxylic_acid(smi)
        print(f"NAME: {name}\nSMILES: {smi}\nResult: {res}\nReason: {reason}\n{'-'*60}")
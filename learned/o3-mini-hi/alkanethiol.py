"""
Classifies: CHEBI:47908 alkanethiol
"""
"""
Classifies: alkanethiol
Definition: An alkanethiol is a compound in which a thiol (-SH) group is attached to an alkyl (i.e. non‐aromatic, “simple”) carbon.
The classifier first rejects molecules that contain extra functional groups (for example, carboxyl, carbonyl, phosphate, sulfonic, ester, amide, or dithiocarbamate groups)
that indicate the –SH function is embedded in a more complex framework.
Then it looks for at least one –SH (thiol) group attached to a non‐aromatic carbon that is not itself directly conjugated to a carbonyl.
"""

from rdkit import Chem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is a simple alkanethiol based on its SMILES string.
    
    A simple alkanethiol is defined as one that contains at least one thiol (-SH) group attached directly to an alkyl (non‐aromatic) carbon and 
    does not contain additional functional groups (e.g. carboxyl, carbonyl (ketone/aldehyde), ester, amide, phosphate, sulfonic, or dithiocarbamate)
    that would complicate the classification.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a simple alkanethiol, False otherwise.
        str: Explanation for the classification decision.
    """
    # Try to parse the SMILES. If invalid, we immediately return an error.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Note: UpdatePropertyCache may be useful.
    mol.UpdatePropertyCache()

    # Disallowed functionality patterns indicating a molecule is more complex than a simple alkanethiol.
    disallowed_smarts = [
        "C(=O)N",         # amide group
        "[CX3](=O)[OX1,H]",# carboxylic acid or ester group
        "[CX3](=O)[#6]",   # ketone or aldehyde (carbonyl adjacent to a carbon)
        "P(=O)(O)O",      # phosphate group
        "S(=O)(=O)[O;H,-]",# sulfonic acid group
        "C(=S)(S)",       # dithiocarbamate motif
    ]
    for smarts in disallowed_smarts:
        pat = Chem.MolFromSmarts(smarts)
        if pat is not None and mol.HasSubstructMatch(pat):
            return False, f"Molecule contains disallowed functionality matching: {smarts}"

    # Now search for valid thiol (-SH) groups.
    # A valid thiol here is defined as:
    #   (1) A sulfur atom (atomic number 16) that has at least one hydrogen (explicit or implicit).
    #   (2) It has exactly one heavy (non-hydrogen) neighbor.
    #   (3) That neighbor is a carbon (atomic number 6) that is not aromatic.
    #   (4) In addition, the carbon neighbor should not be directly linked to a carbonyl (C=O) or an aromatic atom.
    found_valid_thiol = False
    found_thiol = False  # any S with H
    
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16:
            continue  # not sulfur
        # Check that the sulfur carries at least one hydrogen (explicit or implicit)
        if atom.GetTotalNumHs() < 1:
            continue

        # Get heavy (non-hydrogen) neighbors. For a –SH group we expect exactly one.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        if len(heavy_neighbors) != 1:
            continue  # not a simple thiol moiety
        found_thiol = True  # we at least encountered a candidate S

        neighbor = heavy_neighbors[0]
        if neighbor.GetAtomicNum() != 6:
            continue  # sulfur must be bound to a carbon
        if neighbor.GetIsAromatic():
            continue  # disqualify if the carbon is aromatic

        # Check the carbon neighbor's attachments.
        # It should not be directly attached (via a double bond) to an oxygen
        # or be attached to any aromatic atom.
        disqualify = False
        for nbr in neighbor.GetNeighbors():
            if nbr.GetIdx() == atom.GetIdx():
                continue
            bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), nbr.GetIdx())
            if bond is None:
                continue
            # if there is a double bond to an oxygen → carbonyl functionality
            if bond.GetBondTypeAsDouble() == 2 and nbr.GetAtomicNum() == 8:
                disqualify = True
                break
            # if the neighbor is aromatic then do not count this candidate
            if nbr.GetIsAromatic():
                disqualify = True
                break
        if disqualify:
            continue

        # We found a valid –SH group attached to a non-aromatic alkyl carbon.
        found_valid_thiol = True
        break

    if found_valid_thiol:
        return True, "Contains at least one -SH group on a non-aromatic alkyl carbon with no disqualifying adjacent functionalities."
    elif found_thiol:
        # We found a thiol group but its neighboring environment was not acceptable.
        return False, "Found -SH group(s) but attached carbon is conjugated or attached to disallowed functionality."
    else:
        return False, "No valid -SH (thiol) group found in the molecule."

# Example usage when running this script directly.
if __name__ == "__main__":
    test_examples = [
        # True positives (expected to be classified as simple alkanethiols)
        ("SC(C(C)C)C", "3-Methyl-2-butanethiol"),
        ("CC(C)S", "propane-2-thiol"),
        ("SC(C(C)C)CO", "xi-2-Mercapto-3-methyl-1-butanol"),
        ("SC(CCC)C", "2-Pentanethiol"),
        ("SC(C(CO)C)CC", "3-Mercapto-2-methylpentanol"),
        ("SC(CC)CS", "1,2-Butanedithiol"),
        ("SCCCCCCCCCS", "1,9-Nonanedithiol"),
        ("SC(CCCCC)C", "2-Heptanethiol"),
        ("SC(C(S)C)C", "2,3-Butanedithiol"),
        ("CCS", "ethanethiol"),
        ("Cl.SCCN", "Cysteamine hydrochloride"),
        ("SC(CS)(C)C", "2-methylpropane-1,2-dithiol"),
        ("SCCCC", "butanethiol"),
        ("OCCS", "mercaptoethanol"),
        ("SCCCCCCS", "1,6-Hexanedithiol"),
        ("CCSCCS", "2-(ethylsulfanyl)ethanethiol"),
        ("S\\C=C(\\CC)/C", "2-Methyl-1-propenethiol"),
        ("CCCC(C)SCCCO", "4-sulfanyloctan-1-ol"),  # slight variation of the SMILES given
        ("SCC(CC)C", "2-Methyl-1-butanethiol"),
        ("C(C(S)CCO)C", "3-mercaptopentanol"),
        ("SCCCCCCCCS", "1,8-Octanedithiol"),
        ("SCCCCCC", "1-Hexanethiol"),
        ("SCC=C(C)C", "3-Methyl-2-butene-1-thiol"),
        ("SC(CCC)(CO)C", "2-Mercapto-2-methyl-1-pentanol"),
        ("SCC(CCCC)CC", "2-Ethyl-1-hexanethiol"),
        ("SC1CCCC1", "Cyclopentanethiol"),
        # Molecules with extra functionality that should be rejected.
        ("NC(CCS)C(O)=O", "homocysteine"),
        ("S(O)(=O)(=O)CCC(N)CS", "3-Amino-4-sulfanylbutane-1-sulfonic acid"),
        ("SCCC(O)C", "4-Mercapto-2-butanol"),
        ("[O-]C(=O)CC(S)C([O-])=O", "2-mercaptosuccinate"),
        ("[NH3+][C@@H](CCS)C([O-])=O", "L-homocysteine zwitterion"),
        ("SC(CC(=O)C)C", "4-Mercapto-2-pentanone"),
        ("C1CC2C(C1)C3CC2CC3OC(=S)S", "LSM-5150"),
        ("SC(CC)C(=O)C", "3-Mercapto-2-pentanone"),
        ("SC(CCC)CCOC(=O)CCCCC", "3-Mercaptohexyl hexanoate"),
        ("SC(CCCCC(O)=O)CCSCN", "8-[(Aminomethyl)sulfanyl]-6-sulfanyloctanoic acid"),
        ("N[C@@H](CS)C(O)=O", "L-cysteine"),
        ("[O-]C(=O)CS", "thioglycolate(1-)"),
        ("OC(=O)CCCC[C@H](S)CCS", "(S)-dihydrolipoic acid"),
        ("CC(S)C(O)=O", "2-mercaptopropanoic acid"),
        ("SC1(CCC(C(C)C)C(=O)C1)C", "cis- and trans-L-Mercapto-p-menthan-3-one"),
        ("SC(C)C(OCC)=O", "Ethyl 2-mercaptopropionate"),
        ("C([C@H](CS)N)(=O)O.Cl", "L-cysteine hydrochloride"),
        ("S[C@H]1[C@H](O)[C@](O)(C=2C(=O)C=3C(O)=CC=CC3OC2C1)C(=O)OC", "Incarxanthone B"),
        ("SC(C1=CC=CC=C1)C", "1-Phenylethane-1-thiol"),
        ("C(=S)(S)N(CC)CC", "diethyldithiocarbamic acid"),
        ("SCCC(=O)C", "4-Mercapto-2-butanone"),
        ("CC(C)(S)[C@H](N)C(O)=O", "L-penicillamine"),
        ("[NH3+][C@@H](CS)C([O-])=O", "L-cysteine zwitterion"),
        ("[O-]C(=O)CCS", "3-mercaptopropionate"),
        ("C(=S)(S)N(CCCC)CCCC", "dibutyldithiocarbamic acid"),
    ]
    
    for smi, name in test_examples:
        result, explanation = is_alkanethiol(smi)
        print(f"SMILES: {smi:50s} NAME: {name:40s} -> {result} ({explanation})")
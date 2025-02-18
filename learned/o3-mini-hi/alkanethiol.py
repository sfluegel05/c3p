"""
Classifies: CHEBI:47908 alkanethiol
"""
"""
Classifies: alkanethiol
Definition: An alkanethiol is a compound in which a sulfanyl (-SH) group is attached to an alkyl (i.e. non‐aromatic, “simple”) carbon.
This classifier first rejects molecules that contain additional functional groups (e.g. carbonyl, phosphate, sulfonic acid, dithiocarbamate and amide groups)
that indicate the –SH is embedded in a more complex framework. Next it looks for at least one –SH group where the sulfur
(a) carries at least one hydrogen;
(b) is attached to exactly one heavy neighbor which is a carbon;
(c) that carbon is non‐aromatic, sp³ (i.e. “saturated”) and is not directly bound (via double bonds or certain heteroatoms) to groups 
that would indicate additional functionalities (e.g. nitrogen, phosphorus, or a carbonyl bonded oxygen).
"""

from rdkit import Chem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is a simple alkanethiol based on its SMILES string.
    
    A valid alkanethiol contains at least one thiol (-SH) group attached directly to an sp3 (saturated) non‐aromatic alkyl carbon.
    Additionally, the molecule must not contain extra functional groups such as carbonyl, phosphate, sulfonic or dithiocarbamate groups.
    Also, the carbon directly bound to the SH must not be directly attached to disallowed atoms (e.g. nitrogen, phosphorus, halogens)
    or be involved in a C=O bond.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a simple alkanethiol, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    mol.UpdatePropertyCache(strict=False)
    
    # Disallowed functionality patterns.
    # We reject molecules that include groups beyond a simple alkyl framework.
    disallowed_smarts = [
        "[CX3]=[OX1]",      # Carbonyl groups (ketones, aldehydes, carboxylic acids, esters etc.)
        "P(=O)(O)(O)",      # Phosphate/phosphonic acid groups
        "S(=O)(=O)",        # Sulfonic acid groups
        "C(=S)(S)",         # Dithiocarbamate motif
        "C(=O)N",           # Amide groups
    ]
    for smarts in disallowed_smarts:
        pat = Chem.MolFromSmarts(smarts)
        if pat is not None and mol.HasSubstructMatch(pat):
            return False, f"Molecule contains disallowed functionality matching: {smarts}"
    
    # Now search for valid thiol (-SH) candidates.
    valid_candidate = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16:  # must be sulfur
            continue
        # Check that sulfur carries at least one hydrogen (explicit or implicit)
        if atom.GetTotalNumHs() < 1:
            continue
        # Get heavy (non-hydrogen) neighbors
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        if len(heavy_neighbors) != 1:
            continue
        neighbor = heavy_neighbors[0]
        if neighbor.GetAtomicNum() != 6:  # must be carbon
            continue
        if neighbor.GetIsAromatic():    # should be alkyl, not aromatic
            continue
        # Require that the neighboring carbon is sp3 (saturated).
        if neighbor.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        
        # Now inspect the neighbors of the carbon (except the sulfur in question).
        # Allowed atoms for the carbon attached to the thiol: H, C, O, and S.
        # Disallow if the carbon is directly bonded to disallowed atoms (e.g., nitrogen, phosphorus, halogens)
        allowed_atomic_nums = {1, 6, 8, 16}
        disqualify = False
        for nbr in neighbor.GetNeighbors():
            if nbr.GetIdx() == atom.GetIdx():
                continue
            if nbr.GetAtomicNum() not in allowed_atomic_nums:
                disqualify = True
                break
            # Also, if the bond is a double bond to an oxygen then it indicates a carbonyl and must be disqualified.
            bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), nbr.GetIdx())
            if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and nbr.GetAtomicNum() == 8:
                disqualify = True
                break
        if disqualify:
            continue
        
        # Passed all tests: we found a valid simple -SH candidate.
        valid_candidate = True
        break
    
    if valid_candidate:
        return True, "Contains at least one -SH group attached to an sp3 alkyl carbon with no disqualifying adjacent functionalities."
    else:
        return False, "No valid simple alkanethiol group found or molecule contains disallowed functionalities."

# Example usage when running this script directly.
if __name__ == "__main__":
    test_examples = [
        # True positives (should be classified as simple alkanethiols)
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
        ("CCCC(C)SCCCO", "4-sulfanyloctan-1-ol"),
        ("SCC(CC)C", "2-Methyl-1-butanethiol"),
        ("C(C(S)CCO)C", "3-mercaptopentanol"),
        ("SCCCCCCCCS", "1,8-Octanedithiol"),
        ("SCCCCCC", "1-Hexanethiol"),
        ("SCC=C(C)C", "3-Methyl-2-butene-1-thiol"),
        ("SC(CCC)(CO)C", "2-Mercapto-2-methyl-1-pentanol"),
        ("SCC(CCCC)CC", "2-Ethyl-1-hexanethiol"),
        ("SC1CCCC1", "Cyclopentanethiol"),
        # Molecules that should be rejected (for extra functionality or non-alkyl attachment)
        ("SCCC(O)C", "4-Mercapto-2-butanol"),
        ("SC#N", "thiocyanic acid"),
        ("O[C@@H](CS)[C@@H](O)CS", "L-1,4-dithiothreitol"),
        ("SC(C(O)C)C", "3-Mercapto-2-butanol"),
        ("S/C(=N\\C1=C(CC)C=CC=C1CC)/NN", "N1-(2,6-diethylphenyl)hydrazine-1-carbothioamide"),
        ("S(C(S)CCCCC)C", "1-(Methylthio)-1-hexanethiol"),
        ("S1CCC(S)=C1C", "4,5-Dihydro-2-methyl-3-thiophenethiol"),
        ("NCCS", "cysteamine"),
        ("S(C(C(S)C)C)C(C(O)C)C", "3-[(2-Mercapto-1-methylpropyl)thio]-2-butanol"),
        ("SC(C1CCC(=CC1)C)(C)C", "2-(4-Methylcyclohex-3-en-1-yl)propane-2-thiol"),
        ("SC(CC(O)C)(C)C", "(+/-)-4-Mercapto-4-methyl-2-pentanol"),
        ("SC(CCOC=O)(C)C", "3-Mercapto-3-methylbutyl formate"),
        ("O([C@@H]1[C@H]([C@H]([C@@H]([C@H](O1)CO)O)O[C@@H]2[C@H]([C@H]([C@@H]([C@H](O2)CO)O)O)O)C[C@H]3O[C@@H]([C@H]([C@H]([C@@H]3O)O)O)O[C@H]4[C@@H]([C@H]([C@@H](O[C@@H]4CO)OCCCCCCS)N)O", "alpha-D-Man-(1->2)-alpha-D-Man-(1->6)-alpha-D-Man-(1->4)-beta-D-GlcN-O[CH2]6SH"),
        ("S(C(S)CC)CCC", "xi-1-(Propylthio)-1-propanethiol"),
        ("OCC(O)CS", "monothioglycerol"),
        ("C1=NC2=C(C(=N1)N)N=CN2C3[C@H]([C@H](C(O3)CO)O)S", "(3S,4S)-5-(6-aminopurin-9-yl)-2-(hydroxymethyl)-4-mercapto-3-oxolanol"),
        ("SCCCS", "propane-1,3-dithiol"),
        ("OC(CS)C(O)CS", "1,4-dimercaptobutane-2,3-diol"),
        ("CSCCC\\C(S)=N\\O", "4-Methylthiobutylthiohydroximate"),
        ("O[C@H](CS)[C@H](O)CS", "D-1,4-dithiothreitol"),
        ("OCC(S)CS", "dimercaprol"),
        ("S(C(S)CC)C(S)CC", "1,1'-Thiobis-1-propanethiol"),
        ("O[C@@H](CS)[C@H](O)CS", "dithioerythritol"),
        ("SC1CC(OC1C)C", "2,5-Dimethyl-3-mercaptotetrahydrofuran"),
        ("SC(CCOC)(C)C", "4-Methoxy-2-methyl-2-butanethiol"),
    ]
    
    for smi, name in test_examples:
        result, explanation = is_alkanethiol(smi)
        print(f"SMILES: {smi:70s} NAME: {name:40s} -> {result} ({explanation})")
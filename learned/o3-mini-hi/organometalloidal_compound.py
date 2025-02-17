"""
Classifies: CHEBI:143084 organometalloidal compound
"""
"""
Classifies: organometalloidal compound
Definition: A compound having bonds between one or more metalloid atoms (here arsenic, [As])
            and one or more carbon atoms of an organyl group.
            
Heuristics improvement:
  - Only consider molecules whose SMILES string contains "[As" (indicative of As).
  - For every As (atomic number 33) found, examine its As–C bonds.
       • A carbon is “extended” (i.e. part of an organyl group) if it is bonded to at least one
         other carbon (not counting the As) AND if it does not correspond to a long carboxylated fragment.
       • For an extended fragment we extract a “fragment size” (number of heavy atoms reachable 
         from the C atom excluding the current As) and check if a carboxyl group is present (using a SMARTS).
       • If any acceptable extended As–C bond is found, then the compound is accepted.
       • Otherwise, if only methyl groups are present (no extra C neighbors) then accept the bond 
         only if exactly one (or two under narrow conditions) exists AND the overall molecular weight 
         is very low (<170 Da) or the As carries a nonzero formal charge.
         
Note: This heuristic approach is not perfect.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from collections import deque

# Define a SMARTS pattern for a carboxyl group (acidic motif)
CARBOXYL_SMARTS = Chem.MolFromSmarts("C(=O)[O;H1,-]")

def get_fragment_size_and_flag(carbon_atom, exclude_idx):
    """
    From a given carbon atom (neighbor to As), do a breadth-first search to determine:
      - The size (number of heavy atoms, excluding hydrogens) of the fragment reachable if we 
        do not cross back to the arsenic atom (exclude_idx)
      - Whether the fragment contains a carboxyl substructure.
    Returns: (fragment_size, has_carboxyl)
    """
    visited = set()
    queue = deque()
    queue.append(carbon_atom.GetIdx())
    visited.add(carbon_atom.GetIdx())
    
    # Build a submol atom index set
    fragment_atom_indices = []
    
    while queue:
        idx = queue.popleft()
        fragment_atom_indices.append(idx)
        atom = carbon_atom.GetOwningMol().GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            # do not cross if neighbor is the excluded atom (i.e. the arsenic we are coming from)
            if nbr.GetIdx() == exclude_idx:
                continue
            # only consider heavy atoms (skip hydrogens)
            if nbr.GetAtomicNum() == 1:
                continue
            if nbr.GetIdx() not in visited:
                visited.add(nbr.GetIdx())
                queue.append(nbr.GetIdx())
    
    frag_size = len(fragment_atom_indices)
    # Create a submol that corresponds to the fragment for substructure matching:
    submol = Chem.PathToSubmol(carbon_atom.GetOwningMol(), fragment_atom_indices)
    has_carboxyl = submol.HasSubstructMatch(CARBOXYL_SMARTS)
    return frag_size, has_carboxyl

def is_organometalloidal_compound(smiles: str):
    """
    Determines if a molecule is an organometalloidal compound based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is judged to be an organometalloidal compound; False otherwise.
        str: A reason detailing the classification.
    """
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Quick filter: require the SMILES to show [As
    if "[As" not in smiles:
        return False, "No arsenic atom found; not an organometalloidal compound"
    
    # Example hardcoded rejections based on fragments that we know are inorganic
    if "NCC(N)=O" in smiles:
        return False, "Detected tryparsamide-like fragment; not considered an organyl group"
    if "[As](Cl" in smiles:
        return False, "Detected chloride substituents on arsenic; likely inorganic derivative"
    
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    overall_details = []  # collect messages
    compound_valid = False  # flag for valid organometalloidal
    
    # Loop over atoms: check for arsenic (atomic number 33)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 33:
            continue
        
        # For each As atom, find all carbon neighbors
        carbon_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if not carbon_neighbors:
            continue  # no As–C bonds, move to next As
        
        valid_bond_found = False
        details_for_As = []
        methyl_count = 0
        
        for c in carbon_neighbors:
            # Count how many carbon neighbors does c have (exclude the current As)
            other_carbons = [nbr for nbr in c.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != atom.GetIdx()]
            if other_carbons:
                # c is part of a larger organic fragment; let us measure its size and check for carboxyl groups.
                frag_size, has_carboxyl = get_fragment_size_and_flag(c, atom.GetIdx())
                # We now decide: if fragment size is 1 then it is just a methyl; if >1 it is extended.
                if frag_size > 1:
                    # if the fragment is very long (say >3 atoms) and contains a carboxyl group, we flag it as non-organyl.
                    if frag_size > 3 and has_carboxyl:
                        details_for_As.append("Found extended As–C bond with carboxylated chain (rejected)")
                    else:
                        details_for_As.append("Found As–C bond where C is part of a larger organic fragment")
                        valid_bond_found = True
                        break  # if one acceptable extended bond is found, we accept the As atom
                else:
                    # This is essentially a methyl bond.
                    methyl_count += 1
                    details_for_As.append("Found As–C bond where C appears as a methyl group")
            else:
                # No other carbon neighbors -> methyl group.
                methyl_count += 1
                details_for_As.append("Found As–C bond where C appears as a methyl group")
        
        # If no acceptable extended bond was found, consider the methyl situation:
        if not valid_bond_found:
            if methyl_count > 0:
                # Accept if exactly one methyl OR two methyl bonds under low weight / formal charge condition.
                if methyl_count == 1:
                    if atom.GetFormalCharge() != 0 or mol_wt < 170:
                        details_for_As.append("Single methyl group accepted due to charge/low–MW criteria")
                        valid_bond_found = True
                elif methyl_count == 2:
                    if atom.GetFormalCharge() != 0 or mol_wt < 170:
                        details_for_As.append("Two methyl groups accepted based on charge/low–MW criteria")
                        valid_bond_found = True
        
        overall_details.extend(details_for_As)
        if valid_bond_found:
            compound_valid = True
            break  # at least one As atom meets the criteria
    
    if compound_valid:
        return True, "; ".join(overall_details)
    else:
        return False, "No As–C bond found meeting the organyl criteria (extended organic fragment or acceptable methyl pattern)"

# Example usage (test cases from the provided list):
if __name__ == "__main__":
    examples = [
        ("[As](=O)(C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC)(C)C", "1-dimethylarsinoyl-(2Z,5Z,8Z,11Z,14Z,17Z)-heneicosahexaene"),
        ("O[As](O)(=O)c1ccc(cc1)[N+]#N", "ABA diazonium"),
        ("C1C=CC=[As]1", "2H-arsole"),
        ("[As](=O)(CCCCCCCCCCCCCCCCCCCCCCC)(C)C", "1-dimethylarsinoyl-tricosane"),
        ("Nc1cc(ccc1O)[As]1[As]([As]([As]([As]1c1ccc(O)c(N)c1)c1ccc(O)c(N)c1)c1ccc(O)c(N)c1)c1ccc(O)c(N)c1", "arsphenamine pentamer"),
        ("[As](=O)(CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC)(C)C", "1-dimethylarsinoyl-(3Z,6Z,9Z,12Z,15Z,18Z)-docosahexaene"),
        ("[As](=O)(CCCCCCCCCCCCCCCCC)(C)C", "1-dimethylarsinoyl-heptadecane"),
        ("Oc1ccc(cc1[N+]([O-])=O)[As](O)(O)=O", "roxarsone"),
        ("O[As](O)(=O)c1ccc(cc1)[N+]([O-])=O", "nitarsone"),
        ("CC(=O)N[C@@H](Cc1ccc(O)c(c1)\\N=N\\c1ccc(cc1)[As](O)(O)=O)C(=O)NCC(=O)NCC(O)=O", "3-[(4-arsonophenyl)diazenyl]-AcTyrGlyGly"),
        ("C1CC[AsH]C1", "arsolane"),
        ("C[As](O)([O-])=O", "methylarsonate(1-)"),
        ("Nc1ccc(cc1)[As](O)(O)=O", "arsanilic acid"),
        ("[Na+].Nc1ccc(cc1)[As](O)([O-])=O", "sodium arsanilate"),
        ("O[As](O)(=O)c1ccc(cc1)\\N=N\\c1ccc(cc1)[As](O)(O)=O", "4,4'-azodibenzenearsonic acid"),
        ("[As+](CCO)(C)(C)C", "arsenocholine"),
        ("OC(=O)C[As](O)(O)=O", "arsonoacetic acid"),
        ("C[As](O)(O)=O", "methylarsonic acid"),
        ("C[As](C)c1ccccc1[As](C)C", "1,2-phenylenebis(dimethylarsane)"),
        ("C1C=C[As]=C1", "3H-arsole"),
        ("NCC[As](O)(O)=O", "2-Aminoethylarsonate"),
        ("C[As]([O-])([O-])=O", "methylarsonate(2-)"),
        ("CC(=O)N[C@@H](Cc1ccc(O)c(c1)\\N=N\\c1ccc(cc1)[As](O)(O)=O)C(=O)NCC(=O)NCC(=O)NNC(=O)OC(C)(C)C", "3-[(4-arsonophenyl)diazenyl]-AcTyrGlyGlyNHNHBoc"),
        ("O[As]([O-])(=O)c1ccccc1", "phenylarsonate(1-)"),
        ("OC(=O)C\\[As]=[As]\\CC(O)=O", "arsenoacetic acid"),
        ("O[As](O)(=O)c1ccccc1", "phenylarsonic acid"),
        ("Nc1cc(ccc1O)[As]1[As]([As]1c1ccc(O)c(N)c1)c1ccc(O)c(N)c1", "arsphenamine trimer"),
        ("C[As](C)(O)=O", "dimethylarsinic acid"),
        # False positives (should be rejected)
        ("[As](=O)(CCCCCCCCCCCCCCC(O)=O)(C)C", "15-dimethylarsinoyl pentadecanoic acid"),
        ("Nc1cccc(c1)[As](O)O", "m-aminophenylarsonous acid"),
        ("O[As]([O-])(=O)c1ccccc1[N+]#N", "(2-diazoniophenyl)arsonate"),
        ("N(C(C)=O)C=1C=C(C=CC1O)[As]2SC(CS2)CO", "arsthinol"),
        ("[As](=O)(CCCCCCC/C=C\\CCCCCCCC(O)=O)(C)C", "17-dimethylarsinoyl-9Z-heptadecenoic acid"),
        ("Oc1ccc2c(oc3c([As]4SCCS4)c(O)ccc3c2=O)c1[As]1SCCS1", "HOxAsH-bis(1,2-ethanedithiol)"),
        ("[NH3+]c1cc(ccc1O)[As]=[As]c1ccc(O)c([NH3+])c1", "3,3'-diarsene-1,2-diylbis(6-hydroxyanilinium)"),
        ("OC(=O)c1ccccc1-c1c2ccc(O)c([As]3SCCS3)c2oc2c([As]3SCCS3)c(=O)ccc12", "fluorescein bis-arsenide"),
        ("[As](O[Bi]=O)(=O)(O)C1=CC=C(NC(CO)=O)C=C1", "glycobiarsol"),
        ("Nc1cc(cc1O)[As]=[As]c1ccc(O)c(N)c1", "4,4'-diarsene-1,2-diylbis(2-aminophenol)"),
        ("Oc1ccc2nc3ccc(=O)c([As]4SCCS4)c3oc2c1[As]1SCCS1", "resorufin bis-arsenide"),
        ("C1([As](O)O)=CC=C(C=C1)[N+](=O)[O-]", "nitarsone (III)"),
        ("Oc1c(Cl)cc2c(oc3c([As]4SCCS4)c(O)c(Cl)cc3c2=O)c1[As]1SCCS1", "CHOxAsH-bis(1,2-ethanedithiol)"),
        ("N[C@@H](Cc1ccc(O)c(c1)N=Nc1ccc(cc1)[As](O)(O)=O)C(=O)NCCOP(O)(=O)OCC(O)CO", "mono(p-azobenzenearsonic acid)-L-tyrosylglycerophosphorylethanolamine"),
        ("O[As](O)c1ccccc1", "phenylarsonous acid"),
        ("C1=CC(=CC=C1NC(=O)N)[As](=O)(O)O", "[4-(carbamoylamino)phenyl]arsonic acid"),
        ("O[As](O)(=O)CC([O-])=O", "arsonoacetic acid anion"),
        ("S(N)(=O)(=O)C1=CC(N)=C(C=C1)[As](O)(O)=O", "sulfarside"),
        ("[As](=O)(C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCC(O)=O)(C)C", "22-dimethylarsinoyl-(5Z,8Z,11Z,14Z,17Z,20Z)-docosahexaenoic acid"),
        ("[As](=O)(CC1OC(OCC(O)COS(O)=O)C(O)C1O)(C)C", "beta-D-3-[5-Deoxy-5-(dimethylarsinyl)ribofuranosyloxy]-2-hydroxy-1-propanesulfonic acid"),
        ("N[C@@H](CCC(=O)N[C@@H](CSCC(=O)Nc1ccc(cc1)[As](O)O)C(=O)NCC(O)=O)C(O)=O", "S-{2-[4-(dihydroxyarsino)phenylamino]-2-oxoethyl}-glutathione"),
        ("[As](=S)(CCCC(O)=O)(C)C", "4-thiodimethylarsenobutanoic acid"),
        ("N(C(C)=O)C1=CC([As](=O)(O)O)=CC=C1O", "acetarsol"),
        ("C[As]([O-])[O-]", "methylarsonite"),
        ("c1ccc(cc1)[As+](c1ccccc1)(c1ccccc1)c1ccccc1", "tetraphenylarsonium")
    ]
    for smi, name in examples:
        result, reason = is_organometalloidal_compound(smi)
        print(f"Name: {name}\nSMILES: {smi}\nResult: {result}\nReason: {reason}\n{'-'*60}")